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
from ffcx.ir.integral import CommonExpressionIR, TensorPart, compute_integral_ir
from ffcx.ir.representationutils import QuadratureRule, create_quadrature_points_and_weights

logger = logging.getLogger("ffcx")


def basix_cell_from_string(string: str) -> basix.CellType:
    """Convert a string to a Basix CellType."""
    # Note: vertex -> point, rest identity
    match string:
        case "vertex":
            return basix.CellType.point
        case "interval":
            return basix.CellType.interval
        case "triangle":
            return basix.CellType.triangle
        case "tetrahedron":
            return basix.CellType.tetrahedron
        case "quadrilateral":
            return basix.CellType.quadrilateral
        case "hexahedron":
            return basix.CellType.hexahedron
        case "prism":
            return basix.CellType.prism
        case "pyramid":
            return basix.CellType.pyramid
        case _:
            raise KeyError(f"Can not map '{string}' to a basix type.")


class FormIR(typing.NamedTuple):
    """Intermediate representation of a form."""

    id: int
    name: str
    signature: str
    rank: int
    num_coefficients: int
    name_from_uflfile: str
    original_coefficient_positions: list[int]
    coefficient_names: list[str]
    num_constants: int
    constant_ranks: list[int]
    constant_shapes: list[list[int]]
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


class IntegralIR(typing.NamedTuple):
    """Intermediate representation of an integral."""

    expression: CommonExpressionIR
    rank: int
    enabled_coefficients: list[bool]
    part: TensorPart


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
    options: dict[str, npt.DTypeLike | int | float | bool],
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
        form_names[fd_index] = naming.form_name(fd.original_form, fd_index, prefix)  # type: ignore
        for itg_index, itg_data in enumerate(fd.integral_data):  # type: ignore
            integral_names[(fd_index, itg_index)] = naming.integral_name(
                fd.original_form,  # type: ignore
                itg_data.integral_type,
                fd_index,
                itg_data.subdomain_id,
                prefix,
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
    diagonalise = TensorPart.from_str(str(options["part"]))
    ir_forms = [
        _compute_form_ir(
            fd,
            i,
            prefix,
            form_names,
            integral_names,
            integral_domains,
            object_names,
            diagonalise,
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
    form_data,
    form_index,
    element_numbers,
    integral_names,
    options,
    visualise,
) -> list[IntegralIR]:
    """Compute intermediate representation for form integrals."""
    _entity_types = {
        "cell": "cell",
        "exterior_facet": "facet",
        "interior_facet": "facet",
        "vertex": "vertex",
        "custom": "cell",
        "ridge": "ridge",
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
            "part": TensorPart.from_str(options["part"]),
        }
        diagonalise = False
        if form_data.rank == 2 and ir["part"] == TensorPart.diagonal:
            diagonalise = True
            ir["rank"] = 1
            assert form_data.argument_elements[0] == form_data.argument_elements[1], (
                "Can only diagonalise forms with identical arguments."
            )

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

        # Modify output tensor shape if diagonalizing
        if diagonalise:
            expression_ir["tensor_shape"] = expression_ir["tensor_shape"][:1]

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
                if "facet" in integral_type:
                    facet_types = basix.cell.subentity_types(cell_type)[-2]
                    assert len(set(facet_types)) == 1
                    cell_type = facet_types[0]
                elif integral_type == "ridge":
                    ridge_types = basix.cell.subentity_types(cell_type)[-3]
                    assert len(set(ridge_types)) == 1
                    cell_type = ridge_types[0]

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
    form_data,
    form_id,
    prefix,
    form_names,
    integral_names,
    integral_domains,
    object_names,
    tensor_part: TensorPart,
) -> FormIR:
    """Compute intermediate representation of form."""
    logger.info(f"Computing IR for form {form_id}")

    # Store id
    ir = {"id": form_id}

    # Compute common data
    ir["name"] = form_names[form_id]

    ir["signature"] = form_data.original_form.signature()
    args = form_data.original_form.arguments()
    if tensor_part == TensorPart.diagonal and len(args) == 2:
        assert args[0].ufl_function_space() == args[1].ufl_function_space(), (
            "Can only diagonalise forms with identical arguments."
        )
        ir["rank"] = 1
    else:
        ir["rank"] = len(form_data.original_form.arguments())

    ir["num_coefficients"] = len(form_data.reduced_coefficients)

    ir["coefficient_names"] = [
        object_names.get(id(obj), f"w{j}") for j, obj in enumerate(form_data.reduced_coefficients)
    ]

    ir["num_constants"] = len(form_data.original_form.constants())
    ir["constant_ranks"] = [len(obj.ufl_shape) for obj in form_data.original_form.constants()]
    ir["constant_shapes"] = [obj.ufl_shape for obj in form_data.original_form.constants()]

    ir["constant_names"] = [
        object_names.get(id(obj), f"c{j}")
        for j, obj in enumerate(form_data.original_form.constants())
    ]

    ir["original_coefficient_positions"] = form_data.original_coefficient_positions

    ir["finite_element_hashes"] = [
        e.basix_hash() for e in form_data.argument_elements + form_data.coefficient_elements
    ]

    form_name = object_names.get(id(form_data.original_form), form_id)

    ir["name_from_uflfile"] = f"form_{prefix}_{form_name}"

    # Store names of integrals and subdomain_ids for this form, grouped
    # by integral types since form points to all integrals it contains,
    # it has to know their names for codegen phase
    ufcx_integral_types = ("cell", "exterior_facet", "interior_facet", "vertex", "ridge")
    ir["subdomain_ids"] = {itg_type: [] for itg_type in ufcx_integral_types}
    ir["integral_names"] = {itg_type: [] for itg_type in ufcx_integral_types}
    ir["integral_domains"] = {itg_type: [] for itg_type in ufcx_integral_types}
    for itg_index, itg_data in enumerate(form_data.integral_data):
        # UFL is using "otherwise" for default integrals (over whole mesh)
        # but FFCx needs integers, so otherwise = -1
        integral_type = itg_data.integral_type
        subdomain_ids = [sid if sid != "otherwise" else -1 for sid in itg_data.subdomain_id]

        if min(subdomain_ids) < -1:
            raise ValueError("Integral subdomain IDs must be non-negative.")
        ir["subdomain_ids"][integral_type] += subdomain_ids
        for _ in range(len(subdomain_ids)):
            iname = integral_names[(form_id, itg_index)]
            ir["integral_names"][integral_type] += [iname]
            ir["integral_domains"][integral_type] += [integral_domains[iname]]

    return FormIR(**ir)


def _compute_expression_ir(
    expr,
    index,
    prefix,
    analysis,
    options,
    visualise,
    object_names,
):
    """Compute intermediate representation of an Expression."""
    logger.info(f"Computing IR for Expression {index}")

    # Compute representation
    ir = {}
    base_ir = {}
    original_expr = (expr[2], expr[1])

    base_ir["name"] = naming.expression_name(original_expr, prefix)

    original_expr = expr[2]
    points = expr[1]
    expr = expr[0]

    try:
        cell = ufl.domain.extract_unique_domain(expr).ufl_cell()
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
    arguments = ufl.algorithms.extract_arguments(expr)
    argument_elements = tuple(f.ufl_function_space().ufl_element() for f in arguments)
    argument_dimensions = [element_dimensions[element] for element in argument_elements]

    tensor_shape = argument_dimensions
    base_ir["tensor_shape"] = tensor_shape

    base_ir["shape"] = list(expr.ufl_shape)

    coefficients = ufl.algorithms.extract_coefficients(expr)
    coefficient_numbering = {}
    for i, coeff in enumerate(coefficients):
        coefficient_numbering[coeff] = i

    # Add coefficient numbering to IR
    base_ir["coefficient_numbering"] = coefficient_numbering

    original_coefficient_positions = []
    original_coefficients = ufl.algorithms.extract_coefficients(original_expr)
    for coeff in coefficients:
        original_coefficient_positions.append(original_coefficients.index(coeff))

    ir["coefficient_names"] = [
        object_names.get(id(obj), f"w{j}") for j, obj in enumerate(coefficients)
    ]

    ir["constant_names"] = [
        object_names.get(id(obj), f"c{j}")
        for j, obj in enumerate(ufl.algorithms.analysis.extract_constants(expr))
    ]

    expr_name = object_names.get(id(original_expr), index)
    ir["name_from_uflfile"] = f"expression_{prefix}_{expr_name}"

    if len(argument_elements) > 1:
        raise RuntimeError("Expression with more than one Argument not implemented.")

    ir["original_coefficient_positions"] = original_coefficient_positions

    coefficient_elements = tuple(f.ufl_element() for f in coefficients)

    offsets = {}
    _offset = 0
    for i, el in enumerate(coefficient_elements):
        offsets[coefficients[i]] = _offset
        _offset += element_dimensions[el]

    # Copy offsets also into IR
    base_ir["coefficient_offsets"] = offsets

    base_ir["integral_type"] = "expression"
    if cell is not None:
        if (tdim := cell.topological_dimension()) == (pdim := points.shape[1]):
            base_ir["entity_type"] = "cell"
        elif tdim - 1 == pdim:
            base_ir["entity_type"] = "facet"
        else:
            raise ValueError(
                f"Expression on domain with topological dimension {tdim}"
                + f"with points of dimension {pdim} not supported."
            )
    else:
        # For spatially invariant expressions, all expressions are evaluated in the cell
        base_ir["entity_type"] = "cell"

    # Build offsets for Constants
    original_constant_offsets = {}
    _offset = 0
    for constant in ufl.algorithms.analysis.extract_constants(original_expr):
        original_constant_offsets[constant] = _offset
        _offset += np.prod(constant.ufl_shape, dtype=int)

    base_ir["original_constant_offsets"] = original_constant_offsets

    expr_domain = ufl.domain.extract_unique_domain(expr)
    base_ir["coordinate_element_hash"] = (
        expr_domain.ufl_coordinate_element().basix_hash() if expr_domain is not None else 0
    )

    weights = np.array([1.0] * points.shape[0])
    rule = QuadratureRule(points, weights)
    integrands = {"": {rule: expr}}

    if cell is None:
        assert (
            len(ir["original_coefficient_positions"]) == 0
            and len(base_ir["original_constant_offsets"]) == 0
        )

    expression_ir = compute_integral_ir(
        cell,
        base_ir["integral_type"],
        base_ir["entity_type"],
        integrands,
        tensor_shape,
        options,
        visualise,
    )

    base_ir.update(expression_ir)
    ir["expression"] = CommonExpressionIR(**base_ir)
    return ExpressionIR(**ir)
