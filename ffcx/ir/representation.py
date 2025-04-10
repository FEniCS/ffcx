# Copyright (C) 2009-2020 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, Matthew W. Scroggs, Chris Richardson, and others
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 2: Code representation.

Module computes intermediate representations of forms. For each UFC
function, we extract the data needed for code generation at a later
stage.

The representation should conform strictly to the naming and order of
functions in UFC. Thus, for code generation of the function "foo", one
should only need to use the data stored in the intermediate
representation under the key "foo".
"""

from __future__ import annotations

import itertools
import logging
import typing
import warnings

import basix
import numpy as np
import numpy.typing as npt
import ufl
from ufl.classes import Integral
from ufl.sorting import sorted_expr_sum

from ffcx import naming
from ffcx.analysis import UFLData
from ffcx.definitions import entity_types
from ffcx.ir.integral import compute_integral_ir
from ffcx.ir.representationutils import QuadratureRule, create_quadrature_points_and_weights

logger = logging.getLogger("ffcx")


def basix_cell_from_string(string: str) -> basix.CellType:
    """Convert a string to a Basix CellType."""
    # Note: vertex -> point, rest identity
    if string == "vertex":
        return basix.CellType.point
    elif string == "interval":
        return basix.CellType.interval
    elif string == "triangle":
        return basix.CellType.triangle
    elif string == "tetrahedron":
        return basix.CellType.tetrahedron
    elif string == "quadrilateral":
        return basix.CellType.quadrilateral
    elif string == "hexahedron":
        return basix.CellType.hexahedron
    elif string == "prism":
        return basix.CellType.prism
    elif string == "pyramid":
        return basix.CellType.pyramid
    else:
        raise KeyError(f"Can not map '{string}' to a basix type.")

    # TODO: Replace on 31 Oct 2025 with (Python 3.9 does not support match):
    # match string:
    #     case "vertex":
    #         return basix.CellType.point
    #     case "interval":
    #         return basix.CellType.interval
    #     case "triangle":
    #         return basix.CellType.triangle
    #     case "tetrahedron":
    #         return basix.CellType.tetrahedron
    #     case "quadrilateral":
    #         return basix.CellType.quadrilateral
    #     case "hexahedron":
    #         return basix.CellType.hexahedron
    #     case "prism":
    #         return basix.CellType.prism
    #     case "pyramid":
    #         return basix.CellType.pyramid
    #     case _:
    #         raise KeyError(f"Can not map '{string}' to a basix type.")


class FormIR(typing.NamedTuple):
    """Intermediate representation of a form."""

    id: int
    name: str
    signature: str
    rank: int
    num_coefficients: int
    num_constants: int
    name_from_uflfile: str
    original_coefficient_positions: list[int]
    coefficient_names: list[str]
    constant_names: list[str]
    finite_element_hashes: list[int]
    integral_names: dict[str, list[str]]
    integral_domains: dict[str, list[basix.CellType]]
    subdomain_ids: dict[str, list[int]]


class QuadratureIR(typing.NamedTuple):
    """Intermediate representation of a quadrature rule."""

    cell_shape: str
    points: npt.NDArray[np.float64]
    weights: npt.NDArray[np.float64]


class CommonExpressionIR(typing.NamedTuple):
    """Common-ground for IntegralIR and ExpressionIR."""

    integral_type: str
    entity_type: entity_types
    tensor_shape: list[int]
    coefficient_numbering: dict[ufl.Coefficient, int]
    coefficient_offsets: dict[ufl.Coefficient, int]
    original_constant_offsets: dict[ufl.Constant, int]
    unique_tables: dict[str, dict[basix.CellType, npt.NDArray[np.float64]]]
    unique_table_types: dict[basix.CellType, dict[str, str]]
    integrand: dict[tuple[basix.CellType, QuadratureRule], dict]
    name: str
    needs_facet_permutations: bool
    shape: list[int]
    coordinate_element_hash: str


class IntegralIR(typing.NamedTuple):
    """Intermediate representation of an integral."""

    expression: CommonExpressionIR
    rank: int
    enabled_coefficients: list[bool]


class ExpressionIR(typing.NamedTuple):
    """Intermediate representation of a DOLFINx Expression."""

    expression: CommonExpressionIR
    original_coefficient_positions: list[int]
    coefficient_names: list[str]
    constant_names: list[str]
    name_from_uflfile: str


class DataIR(typing.NamedTuple):
    """Intermediate representation of data."""

    integrals: list[IntegralIR]
    forms: list[FormIR]
    expressions: list[ExpressionIR]


def compute_ir(
    analysis: UFLData,
    object_names: dict[int, str],
    prefix: str,
    options: dict[str, typing.Union[npt.DTypeLike, int, float]],
    visualise: bool,
) -> DataIR:
    """Compute intermediate representation."""
    logger.info(79 * "*")
    logger.info("Compiler stage 2: Computing intermediate representation of objects")
    logger.info(79 * "*")

    # Compute object names
    # NOTE: This is done here for performance reasons, because repeated
    # calls within each IR computation would be expensive due to UFL
    # signature computations.
    integral_names = {}
    form_names = {}
    for fd_index, fd in enumerate(analysis.form_data):
        form_names[fd_index] = naming.form_name(fd.original_form, fd_index, prefix)
        for itg_index, itg_data in enumerate(fd.integral_data):
            integral_names[(fd_index, itg_index)] = naming.integral_name(
                fd.original_form, itg_data.integral_type, fd_index, itg_data.subdomain_id, prefix
            )

    irs = [
        _compute_integral_ir(
            fd,
            i,
            analysis.element_numbers,
            integral_names,
            options,
            visualise,
        )
        for (i, fd) in enumerate(analysis.form_data)
    ]
    ir_integrals = list(itertools.chain(*irs))

    integral_domains = {
        i.expression.name: set(j[0] for j in i.expression.integrand.keys()) for a in irs for i in a
    }

    ir_forms = [
        _compute_form_ir(
            fd,
            i,
            prefix,
            form_names,
            integral_names,
            integral_domains,
            object_names,
        )
        for (i, fd) in enumerate(analysis.form_data)
    ]

    ir_expressions = [
        _compute_expression_ir(
            expr,
            i,
            prefix,
            analysis,
            options,
            visualise,
            object_names,
        )
        for i, expr in enumerate(analysis.expressions)
    ]

    return DataIR(
        integrals=ir_integrals,
        forms=ir_forms,
        expressions=ir_expressions,
    )


def _compute_integral_ir(
    form_data: ufl.classes.Form,
    form_index: int,
    element_numbers: dict[typing.Any, int],
    integral_names: dict[tuple[int, int], str],
    options: dict[str, typing.Union[npt.DTypeLike, int, float, bool]],
    visualise: bool,
) -> list[IntegralIR]:
    """Compute intermediate representation for form integrals."""
    _entity_types = {
        "cell": "cell",
        "exterior_facet": "facet",
        "interior_facet": "facet",
        "vertex": "vertex",
        "custom": "cell",
    }

    # Iterate over groups of integrals
    irs = []
    for itg_data_index, itg_data in enumerate(form_data.integral_data):
        logger.info(f"Computing IR for integral in integral group {itg_data_index}")
        expression_ir = {}

        # Compute representation
        entity_type = _entity_types[itg_data.integral_type]
        ufl_cell = itg_data.domain.ufl_cell()
        cell_type = basix_cell_from_string(ufl_cell.cellname())
        tdim = ufl_cell.topological_dimension()
        assert all(tdim == itg.ufl_domain().topological_dimension() for itg in itg_data.integrals)

        expression_ir = {
            "integral_type": itg_data.integral_type,
            "entity_type": entity_type,
            "shape": (),
            "coordinate_element_hash": itg_data.domain.ufl_coordinate_element().basix_hash(),
        }
        ir = {
            "rank": form_data.rank,
            "enabled_coefficients": itg_data.enabled_coefficients,
        }

        # Get element space dimensions
        unique_elements = element_numbers.keys()
        element_dimensions = {
            element: element.dim + element.num_global_support_dofs for element in unique_elements
        }

        # Create dimensions of primary indices, needed to reset the argument
        # 'A' given to tabulate_tensor() by the assembler.
        argument_dimensions = [
            element_dimensions[element] for element in form_data.argument_elements
        ]

        # Compute shape of element tensor
        if expression_ir["integral_type"] == "interior_facet":
            expression_ir["tensor_shape"] = [2 * dim for dim in argument_dimensions]
        else:
            expression_ir["tensor_shape"] = argument_dimensions

        integral_type = itg_data.integral_type

        # Group integrands with the same quadrature rule
        grouped_integrands: dict[
            basix.CellType, dict[QuadratureRule, list[ufl.core.expr.Expr]]
        ] = {}
        use_sum_factorization = options["sum_factorization"] and itg_data.integral_type == "cell"
        for integral in itg_data.integrals:
            md = integral.metadata() or {}
            scheme = md["quadrature_rule"]
            tensor_factors = None
            rules = {}
            if scheme == "custom":
                points = md["quadrature_points"]
                weights = md["quadrature_weights"]
                rules[cell_type] = (points, weights, None)
            elif scheme == "vertex":
                # The vertex scheme, i.e., averaging the function value in the
                # vertices and multiplying with the simplex volume, is only of
                # order 1 and inferior to other generic schemes in terms of
                # error reduction. Equation systems generated with the vertex
                # scheme have some properties that other schemes lack, e.g., the
                # mass matrix is a simple diagonal matrix. This may be
                # prescribed in certain cases.

                degree = md["quadrature_degree"]
                if integral_type != "cell":
                    facet_types: list[basix.CellType] = basix.cell.subentity_types(cell_type)[-2]
                    assert len(set(facet_types)) == 1
                    cell_type = facet_types[0]
                if degree > 1:
                    warnings.warn(
                        "Explicitly selected vertex quadrature (degree 1), "
                        f"but requested degree is {degree}."
                    )
                points = basix.cell.geometry(cell_type)
                cell_volume = basix.cell.volume(cell_type)
                weights = np.full(
                    points.shape[0], cell_volume / points.shape[0], dtype=points.dtype
                )
                rules[cell_type] = (points, weights, None)
            else:
                degree = md["quadrature_degree"]
                points, weights, tensor_factors = create_quadrature_points_and_weights(
                    integral_type,
                    ufl_cell,
                    degree,
                    scheme,
                    form_data.argument_elements,
                    use_sum_factorization,
                )
                rules = {
                    basix_cell_from_string(i): (
                        points[i],
                        weights[i],
                        tensor_factors[i] if i in tensor_factors else None,
                    )
                    for i in points
                }

            for cell_type, (points, weights, tensor_factors) in rules.items():
                points = np.asarray(points)
                weights = np.asarray(weights)
                rule = QuadratureRule(points, weights, tensor_factors)

                if cell_type not in grouped_integrands:
                    grouped_integrands[cell_type] = {}
                if rule not in grouped_integrands:
                    grouped_integrands[cell_type][rule] = []
                grouped_integrands[cell_type][rule].append(integral.integrand())
        sorted_integrals: dict[basix.CellType, dict[QuadratureRule, Integral]] = {
            cell_type: {} for cell_type in grouped_integrands
        }
        for cell_type, integrands_by_cell in grouped_integrands.items():
            for rule, integrands in integrands_by_cell.items():
                integrands_summed = sorted_expr_sum(integrands)

                integral_new = Integral(
                    integrands_summed,
                    itg_data.integral_type,
                    itg_data.domain,
                    itg_data.subdomain_id,
                    {},
                    None,
                )
                sorted_integrals[cell_type][rule] = integral_new

        # TODO: See if coefficient_numbering can be removed
        # Build coefficient numbering for UFC interface here, to avoid
        # renumbering in UFL and application of replace mapping
        coefficient_numbering = {}
        for i, f in enumerate(form_data.reduced_coefficients):
            coefficient_numbering[f] = i

        # Add coefficient numbering to IR
        expression_ir["coefficient_numbering"] = coefficient_numbering

        index_to_coeff = sorted([(v, k) for k, v in coefficient_numbering.items()])
        offsets = {}
        width = 2 if integral_type in ("interior_facet") else 1
        _offset = 0
        for k, el in zip(index_to_coeff, form_data.coefficient_elements):
            offsets[k[1]] = _offset
            _offset += width * element_dimensions[el]

        # Copy offsets also into IR
        expression_ir["coefficient_offsets"] = offsets

        # Build offsets for Constants
        original_constant_offsets = {}
        _offset = 0
        for constant in form_data.original_form.constants():
            original_constant_offsets[constant] = _offset
            _offset += np.prod(constant.ufl_shape, dtype=int)

        expression_ir["original_constant_offsets"] = original_constant_offsets

        # Create map from number of quadrature points -> integrand
        integrand_map: dict[basix.CellType, dict[QuadratureRule, ufl.core.expr.Expr]] = {
            cell_type: {rule: integral.integrand() for rule, integral in cell_integrals.items()}
            for cell_type, cell_integrals in sorted_integrals.items()
        }

        # Build more specific intermediate representation
        integral_ir = compute_integral_ir(
            itg_data.domain.ufl_cell(),
            itg_data.integral_type,
            expression_ir["entity_type"],
            integrand_map,
            expression_ir["tensor_shape"],
            options,
            visualise,
        )

        expression_ir.update(integral_ir)

        # Fetch name
        expression_ir["name"] = integral_names[(form_index, itg_data_index)]
        ir["expression"] = CommonExpressionIR(**expression_ir)
        irs.append(IntegralIR(**ir))

    return irs


def _compute_form_ir(
    form_data: ufl.classes.Form,
    form_id: int,
    prefix: str,
    form_names: dict[int, str],
    integral_names: dict[tuple[int, int], str],
    integral_domains: dict[str, set[basix.CellType]],
    object_names: dict[int, str],
) -> FormIR:
    """Compute intermediate representation of form."""
    logger.info(f"Computing IR for form {form_id}")

    # Create dictionary for form IR data
    form_id_value = form_id
    name = form_names[form_id]
    signature = form_data.original_form.signature()
    rank = len(form_data.original_form.arguments())
    num_coefficients = len(form_data.reduced_coefficients)
    num_constants = len(form_data.original_form.constants())

    coefficient_names = [
        object_names.get(id(obj), f"w{j}") for j, obj in enumerate(form_data.reduced_coefficients)
    ]

    constant_names = [
        object_names.get(id(obj), f"c{j}")
        for j, obj in enumerate(form_data.original_form.constants())
    ]

    original_coefficient_positions = form_data.original_coefficient_positions

    finite_element_hashes = [
        e.basix_hash() for e in form_data.argument_elements + form_data.coefficient_elements
    ]

    form_name = object_names.get(id(form_data.original_form), form_id)
    name_from_uflfile = f"form_{prefix}_{form_name}"

    # Store names of integrals and subdomain_ids for this form, grouped
    # by integral types since form points to all integrals it contains,
    # it has to know their names for codegen phase
    ufcx_integral_types = ("cell", "exterior_facet", "interior_facet")
    subdomain_ids: dict[str, list[int]] = {itg_type: [] for itg_type in ufcx_integral_types}
    integral_names_dict: dict[str, list[str]] = {itg_type: [] for itg_type in ufcx_integral_types}
    integral_domains_dict: dict[str, list[basix.CellType]] = {itg_type: [] for itg_type in ufcx_integral_types}
    
    for itg_index, itg_data in enumerate(form_data.integral_data):
        # UFL is using "otherwise" for default integrals (over whole mesh)
        # but FFCx needs integers, so otherwise = -1
        integral_type = itg_data.integral_type
        itg_subdomain_ids = [sid if sid != "otherwise" else -1 for sid in itg_data.subdomain_id]

        if min(itg_subdomain_ids) < -1:
            raise ValueError("Integral subdomain IDs must be non-negative.")
        subdomain_ids[integral_type] += itg_subdomain_ids
        for _ in range(len(itg_subdomain_ids)):
            iname = integral_names[(form_id, itg_index)]
            integral_names_dict[integral_type] += [iname]
            # Convert set to a single element for the list
            domain_set = integral_domains[iname]
            if len(domain_set) == 1:
                domain_element = next(iter(domain_set))
            else:
                # If there are multiple cell types, just pick the first one
                # This is a simplification that may need to be addressed in the future
                domain_element = next(iter(domain_set))
            integral_domains_dict[integral_type] += [domain_element]

    return FormIR(
        id=form_id_value,
        name=name,
        signature=signature,
        rank=rank,
        num_coefficients=num_coefficients,
        num_constants=num_constants,
        name_from_uflfile=name_from_uflfile,
        original_coefficient_positions=original_coefficient_positions,
        coefficient_names=coefficient_names,
        constant_names=constant_names,
        finite_element_hashes=finite_element_hashes,
        integral_names=integral_names_dict,
        integral_domains=integral_domains_dict,
        subdomain_ids=subdomain_ids
    )


def _compute_expression_ir(
    expr: tuple[ufl.core.expr.Expr, npt.NDArray[np.float64], ufl.core.expr.Expr],
    index: int,
    prefix: str,
    analysis: UFLData,
    options: dict[str, typing.Union[npt.DTypeLike, int, float, bool]],
    visualise: bool,
    object_names: dict[int, str],
) -> ExpressionIR:
    """Compute intermediate representation of an Expression."""
    logger.info(f"Computing IR for Expression {index}")

    # Compute representation
    original_expr = (expr[2], expr[1])
    name = naming.expression_name(original_expr, prefix)

    original_expr_val = expr[2]
    points = expr[1]
    expr_val = expr[0]

    try:
        cell = ufl.domain.extract_unique_domain(expr_val).ufl_cell()
    except AttributeError:
        # This case corresponds to a spatially constant expression
        # without any dependencies
        cell = None

    # Prepare dimensions of all unique element in expression, including
    # elements for arguments, coefficients and coordinate mappings
    element_dimensions = {
        element: element.dim + element.num_global_support_dofs
        for element in analysis.unique_elements
    }

    # Extract dimensions for elements of arguments only
    arguments = ufl.algorithms.extract_arguments(expr_val)
    argument_elements = tuple(f.ufl_function_space().ufl_element() for f in arguments)
    argument_dimensions = [element_dimensions[element] for element in argument_elements]

    tensor_shape = argument_dimensions
    shape = list(expr_val.ufl_shape)

    coefficients = ufl.algorithms.extract_coefficients(expr_val)
    coefficient_numbering = {}
    for i, coeff in enumerate(coefficients):
        coefficient_numbering[coeff] = i

    original_coefficient_positions = []
    original_coefficients = ufl.algorithms.extract_coefficients(original_expr_val)
    for coeff in coefficients:
        original_coefficient_positions.append(original_coefficients.index(coeff))

    coefficient_names = [
        object_names.get(id(obj), f"w{j}") for j, obj in enumerate(coefficients)
    ]

    constant_names = [
        object_names.get(id(obj), f"c{j}")
        for j, obj in enumerate(ufl.algorithms.analysis.extract_constants(expr_val))
    ]

    expr_name = object_names.get(id(original_expr_val), index)
    name_from_uflfile = f"expression_{prefix}_{expr_name}"

    if len(argument_elements) > 1:
        raise RuntimeError("Expression with more than one Argument not implemented.")

    coefficient_elements = tuple(f.ufl_element() for f in coefficients)

    offsets = {}
    _offset = 0
    for i, el in enumerate(coefficient_elements):
        offsets[coefficients[i]] = _offset
        _offset += element_dimensions[el]

    integral_type = "expression"
    entity_type: entity_types
    if cell is not None:
        if (tdim := cell.topological_dimension()) == (pdim := points.shape[1]):
            entity_type = "cell"
        elif tdim - 1 == pdim:
            entity_type = "facet"
        else:
            raise ValueError(
                f"Expression on domain with topological dimension {tdim}"
                + f"with points of dimension {pdim} not supported."
            )
    else:
        # For spatially invariant expressions, all expressions are evaluated in the cell
        entity_type = "cell"

    # Build offsets for Constants
    original_constant_offsets = {}
    _offset = 0
    for constant in ufl.algorithms.analysis.extract_constants(original_expr_val):
        original_constant_offsets[constant] = _offset
        _offset += np.prod(constant.ufl_shape, dtype=int)

    coordinate_element_hash = ""
    if cell is not None:
        coordinate_element_hash = (
            ufl.domain.extract_unique_domain(expr_val).ufl_coordinate_element().basix_hash()
        )

    weights = np.array([1.0] * points.shape[0])
    rule = QuadratureRule(points, weights)
    integrands = {"": {rule: expr_val}}

    if cell is None:
        assert (
            len(original_coefficient_positions) == 0
            and len(original_constant_offsets) == 0
        )

    # For expressions, ensure entity_type is one of the allowed values from the Literal type
    # This cast is necessary because mypy doesn't recognize that our assignment ensures the correct type
    entity_type_arg: entity_types = entity_type  # Explicit cast to help mypy

    expression_ir = compute_integral_ir(
        cell,
        integral_type,
        entity_type_arg,
        integrands,
        tensor_shape,
        options,
        visualise,
    )

    # Update with computed data
    common_expr_ir = CommonExpressionIR(
        integral_type=integral_type,
        entity_type=entity_type,
        tensor_shape=tensor_shape,
        coefficient_numbering=coefficient_numbering,
        coefficient_offsets=offsets,
        original_constant_offsets=original_constant_offsets,
        unique_tables=expression_ir["unique_tables"],
        unique_table_types=expression_ir["unique_table_types"],
        integrand=expression_ir["integrand"],
        name=name,
        needs_facet_permutations=expression_ir["needs_facet_permutations"],
        shape=shape,
        coordinate_element_hash=coordinate_element_hash
    )

    return ExpressionIR(
        expression=common_expr_ir,
        original_coefficient_positions=original_coefficient_positions,
        coefficient_names=coefficient_names,
        constant_names=constant_names,
        name_from_uflfile=name_from_uflfile
    )
