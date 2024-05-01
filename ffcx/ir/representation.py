# Copyright (C) 2009-2020 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, Matthew W. Scroggs, Chris Richardson, and others
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 2: Code representation.

Module computes intermediate representations of forms. For each UFC
function, we extract the data needed for code generation at a later stage.

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
import basix.ufl
import numpy as np
import numpy.typing as npt
import ufl
from ufl.classes import Integral
from ufl.sorting import sorted_expr_sum

from ffcx import naming
from ffcx.analysis import UFLData
from ffcx.ir.integral import compute_integral_ir
from ffcx.ir.representationutils import QuadratureRule, create_quadrature_points_and_weights

logger = logging.getLogger("ffcx")


class FormIR(typing.NamedTuple):
    """Intermediate representation of a form."""

    id: int
    name: str
    signature: str
    rank: int
    num_coefficients: int
    num_constants: int
    name_from_uflfile: str
    function_spaces: dict[
        str, tuple[str, str, str, int, basix.CellType, basix.LagrangeVariant, tuple[int]]
    ]
    original_coefficient_position: list[int]
    coefficient_names: list[str]
    constant_names: list[str]
    finite_element_hashes: list[int]
    integral_names: dict[str, list[str]]
    subdomain_ids: dict[str, list[int]]


class QuadratureIR(typing.NamedTuple):
    """Intermediate representation of a quadrature rule."""

    cell_shape: str
    points: npt.NDArray[np.float64]
    weights: npt.NDArray[np.float64]


class IntegralIR(typing.NamedTuple):
    """Intermediate representation of an integral."""

    integral_type: str
    rank: int
    entitytype: str
    enabled_coefficients: list[bool]
    tensor_shape: list[int]
    coefficient_numbering: dict[ufl.Coefficient, int]
    coefficient_offsets: dict[ufl.Coefficient, int]
    original_constant_offsets: dict[ufl.Constant, int]
    unique_tables: dict[str, npt.NDArray[np.float64]]
    unique_table_types: dict[str, str]
    integrand: dict[QuadratureRule, dict]
    name: str
    needs_facet_permutations: bool
    coordinate_element_hash: int


class ExpressionIR(typing.NamedTuple):
    """Intermediate representation of an expression."""

    name: str
    options: dict
    unique_tables: dict[str, npt.NDArray[np.float64]]
    unique_table_types: dict[str, str]
    integrand: dict[QuadratureRule, dict]
    coefficient_numbering: dict[ufl.Coefficient, int]
    coefficient_offsets: dict[ufl.Coefficient, int]
    integral_type: str
    entitytype: str
    tensor_shape: list[int]
    expression_shape: list[int]
    original_constant_offsets: dict[ufl.Constant, int]
    points: npt.NDArray[np.float64]
    coefficient_names: list[str]
    constant_names: list[str]
    needs_facet_permutations: bool
    function_spaces: dict[
        str, tuple[str, str, str, int, basix.CellType, basix.LagrangeVariant, tuple[int]]
    ]
    name_from_uflfile: str
    original_coefficient_positions: list[int]


class DataIR(typing.NamedTuple):
    """Intermediate representation of data."""

    integrals: list[IntegralIR]
    forms: list[FormIR]
    expressions: list[ExpressionIR]


def compute_ir(
    analysis: UFLData,
    object_names: dict[int, str],
    prefix: str,
    options: dict[str, npt.DTypeLike | int | float],
    visualise: bool,
) -> DataIR:
    """Compute intermediate representation."""
    logger.info(79 * "*")
    logger.info("Compiler stage 2: Computing intermediate representation of objects")
    logger.info(79 * "*")

    # Compute object names
    # NOTE: This is done here for performance reasons, because repeated calls
    # within each IR computation would be expensive due to UFL signature computations
    finite_element_hashes = {e: e.basix_hash() for e in analysis.unique_elements}
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
            finite_element_hashes,
            options,
            visualise,
        )
        for (i, fd) in enumerate(analysis.form_data)
    ]
    ir_integrals = list(itertools.chain(*irs))

    ir_forms = [
        _compute_form_ir(
            fd,
            i,
            prefix,
            form_names,
            integral_names,
            analysis.element_numbers,
            finite_element_hashes,
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
            finite_element_hashes,
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
    finite_element_hashes,
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
    }

    # Iterate over groups of integrals
    irs = []
    for itg_data_index, itg_data in enumerate(form_data.integral_data):
        logger.info(f"Computing IR for integral in integral group {itg_data_index}")

        # Compute representation
        entitytype = _entity_types[itg_data.integral_type]
        cell = itg_data.domain.ufl_cell()
        cellname = cell.cellname()
        tdim = cell.topological_dimension()
        assert all(tdim == itg.ufl_domain().topological_dimension() for itg in itg_data.integrals)

        ir = {
            "integral_type": itg_data.integral_type,
            "rank": form_data.rank,
            "entitytype": entitytype,
            "enabled_coefficients": itg_data.enabled_coefficients,
            "coordinate_element_hash": finite_element_hashes[
                itg_data.domain.ufl_coordinate_element()
            ],
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
        if ir["integral_type"] == "interior_facet":
            ir["tensor_shape"] = [2 * dim for dim in argument_dimensions]
        else:
            ir["tensor_shape"] = argument_dimensions

        integral_type = itg_data.integral_type
        cell = itg_data.domain.ufl_cell()

        # Group integrands with the same quadrature rule
        grouped_integrands: dict[QuadratureRule, list[ufl.core.expr.Expr]] = {}
        use_sum_factorization = options["sum_factorization"] and itg_data.integral_type == "cell"
        for integral in itg_data.integrals:
            md = integral.metadata() or {}
            scheme = md["quadrature_rule"]
            tensor_factors = None
            if scheme == "custom":
                points = md["quadrature_points"]
                weights = md["quadrature_weights"]
            elif scheme == "vertex":
                # FIXME: Could this come from basix?

                # The vertex scheme, i.e., averaging the function value in the
                # vertices and multiplying with the simplex volume, is only of
                # order 1 and inferior to other generic schemes in terms of
                # error reduction. Equation systems generated with the vertex
                # scheme have some properties that other schemes lack, e.g., the
                # mass matrix is a simple diagonal matrix. This may be
                # prescribed in certain cases.

                degree = md["quadrature_degree"]
                if integral_type != "cell":
                    facet_types = cell.facet_types()
                    assert len(facet_types) == 1
                    cellname = facet_types[0].cellname()
                if degree > 1:
                    warnings.warn(
                        "Explicitly selected vertex quadrature (degree 1), "
                        f"but requested degree is {degree}."
                    )
                if cellname == "tetrahedron":
                    points, weights = (
                        np.array(
                            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
                        ),
                        np.array([1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0]),
                    )
                elif cellname == "triangle":
                    points, weights = (
                        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                        np.array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]),
                    )
                elif cellname == "interval":
                    # Trapezoidal rule
                    points, weights = (np.array([[0.0], [1.0]]), np.array([1.0 / 2.0, 1.0 / 2.0]))
                elif cellname == "quadrilateral":
                    points, weights = (
                        np.array([[0.0, 0], [1.0, 0.0], [0.0, 1.0], [1.0, 1]]),
                        np.array([1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0]),
                    )
                elif cellname == "hexahedron":
                    points, weights = (
                        np.array(
                            [
                                [0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0],
                                [0.0, 1.0, 0.0],
                                [1.0, 1.0, 0.0],
                                [0.0, 0.0, 1.0],
                                [1.0, 0.0, 1.0],
                                [0.0, 1.0, 1.0],
                                [1.0, 1.0, 1.0],
                            ]
                        ),
                        np.array(
                            [
                                1.0 / 8.0,
                                1.0 / 8.0,
                                1.0 / 8.0,
                                1.0 / 8.0,
                                1.0 / 8.0,
                                1.0 / 8.0,
                                1.0 / 8.0,
                                1.0 / 8.0,
                            ]
                        ),
                    )
                else:
                    raise RuntimeError(f"Vertex scheme is not supported for cell: {cellname}")
            else:
                degree = md["quadrature_degree"]
                points, weights, tensor_factors = create_quadrature_points_and_weights(
                    integral_type,
                    cell,
                    degree,
                    scheme,
                    form_data.argument_elements,
                    use_sum_factorization,
                )

            points = np.asarray(points)
            weights = np.asarray(weights)
            rule = QuadratureRule(points, weights, tensor_factors)

            if rule not in grouped_integrands:
                grouped_integrands[rule] = []
            grouped_integrands[rule].append(integral.integrand())
        sorted_integrals: dict[QuadratureRule, Integral] = {}
        for rule, integrands in grouped_integrands.items():
            integrands_summed = sorted_expr_sum(integrands)

            integral_new = Integral(
                integrands_summed,
                itg_data.integral_type,
                itg_data.domain,
                itg_data.subdomain_id,
                {},
                None,
            )
            sorted_integrals[rule] = integral_new

        # TODO: See if coefficient_numbering can be removed
        # Build coefficient numbering for UFC interface here, to avoid
        # renumbering in UFL and application of replace mapping
        coefficient_numbering = {}
        for i, f in enumerate(form_data.reduced_coefficients):
            coefficient_numbering[f] = i

        # Add coefficient numbering to IR
        ir["coefficient_numbering"] = coefficient_numbering

        index_to_coeff = sorted([(v, k) for k, v in coefficient_numbering.items()])
        offsets = {}
        width = 2 if integral_type in ("interior_facet") else 1
        _offset = 0
        for k, el in zip(index_to_coeff, form_data.coefficient_elements):
            offsets[k[1]] = _offset
            _offset += width * element_dimensions[el]

        # Copy offsets also into IR
        ir["coefficient_offsets"] = offsets

        # Build offsets for Constants
        original_constant_offsets = {}
        _offset = 0
        for constant in form_data.original_form.constants():
            original_constant_offsets[constant] = _offset
            _offset += np.prod(constant.ufl_shape, dtype=int)

        ir["original_constant_offsets"] = original_constant_offsets

        # Create map from number of quadrature points -> integrand
        integrand_map: dict[QuadratureRule, ufl.core.expr.Expr] = {
            rule: integral.integrand() for rule, integral in sorted_integrals.items()
        }

        # Build more specific intermediate representation
        integral_ir = compute_integral_ir(
            itg_data.domain.ufl_cell(),
            itg_data.integral_type,
            ir["entitytype"],
            integrand_map,
            ir["tensor_shape"],
            options,
            visualise,
        )

        ir.update(integral_ir)

        # Fetch name
        ir["name"] = integral_names[(form_index, itg_data_index)]

        irs.append(IntegralIR(**ir))

    return irs


def _compute_form_ir(
    form_data,
    form_id,
    prefix,
    form_names,
    integral_names,
    element_numbers,
    finite_element_hashes,
    object_names,
) -> FormIR:
    """Compute intermediate representation of form."""
    logger.info(f"Computing IR for form {form_id}")

    # Store id
    ir = {"id": form_id}

    # Compute common data
    ir["name"] = form_names[form_id]

    ir["signature"] = form_data.original_form.signature()

    ir["rank"] = len(form_data.original_form.arguments())
    ir["num_coefficients"] = len(form_data.reduced_coefficients)
    ir["num_constants"] = len(form_data.original_form.constants())

    ir["coefficient_names"] = [
        object_names.get(id(obj), f"w{j}") for j, obj in enumerate(form_data.reduced_coefficients)
    ]

    ir["constant_names"] = [
        object_names.get(id(obj), f"c{j}")
        for j, obj in enumerate(form_data.original_form.constants())
    ]

    ir["original_coefficient_position"] = form_data.original_coefficient_positions

    ir["finite_element_hashes"] = [
        finite_element_hashes[e]
        for e in form_data.argument_elements + form_data.coefficient_elements
    ]

    fs = {}
    for function in form_data.original_form.arguments() + tuple(form_data.reduced_coefficients):
        name = object_names.get(id(function), str(function))
        if not str(name).isidentifier():
            raise ValueError(f'Function name "{name}" must be a valid object identifier.')
        el = function.ufl_function_space().ufl_element()
        space = function.ufl_function_space()
        domain = space.ufl_domain()
        cmap = domain.ufl_coordinate_element()
        # Default point spacing for CoordinateElement is equispaced
        if not isinstance(cmap, basix.ufl._ElementBase) and cmap.variant() is None:
            cmap._sub_element._variant = "equispaced"
        family = cmap.family_name
        degree = cmap.degree
        value_shape = space.value_shape
        fs[name] = (
            finite_element_hashes[el],
            family,
            degree,
            cmap.cell_type,
            cmap.lagrange_variant,
            value_shape,
        )

    form_name = object_names.get(id(form_data.original_form), form_id)

    ir["function_spaces"] = fs
    ir["name_from_uflfile"] = f"form_{prefix}_{form_name}"

    # Store names of integrals and subdomain_ids for this form, grouped
    # by integral types since form points to all integrals it contains,
    # it has to know their names for codegen phase
    ir["integral_names"] = {}
    ir["subdomain_ids"] = {}
    ufcx_integral_types = ("cell", "exterior_facet", "interior_facet")
    ir["subdomain_ids"] = {itg_type: [] for itg_type in ufcx_integral_types}
    ir["integral_names"] = {itg_type: [] for itg_type in ufcx_integral_types}
    for itg_index, itg_data in enumerate(form_data.integral_data):
        # UFL is using "otherwise" for default integrals (over whole mesh)
        # but FFCx needs integers, so otherwise = -1
        integral_type = itg_data.integral_type
        subdomain_ids = [sid if sid != "otherwise" else -1 for sid in itg_data.subdomain_id]

        if min(subdomain_ids) < -1:
            raise ValueError("Integral subdomain IDs must be non-negative.")
        ir["subdomain_ids"][integral_type] += subdomain_ids
        for _ in range(len(subdomain_ids)):
            ir["integral_names"][integral_type] += [integral_names[(form_id, itg_index)]]

    return FormIR(**ir)


def _compute_expression_ir(
    expression,
    index,
    prefix,
    analysis,
    options,
    visualise,
    object_names,
    finite_element_hashes,
):
    """Compute intermediate representation of expression."""
    logger.info(f"Computing IR for expression {index}")

    # Compute representation
    ir = {}

    original_expression = (expression[2], expression[1])

    ir["name"] = naming.expression_name(original_expression, prefix)

    original_expression = expression[2]
    points = expression[1]
    expression = expression[0]

    try:
        cell = ufl.domain.extract_unique_domain(expression).ufl_cell()
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
    arguments = ufl.algorithms.extract_arguments(expression)
    argument_elements = tuple(f.ufl_function_space().ufl_element() for f in arguments)
    argument_dimensions = [element_dimensions[element] for element in argument_elements]

    tensor_shape = argument_dimensions
    ir["tensor_shape"] = tensor_shape

    ir["expression_shape"] = list(expression.ufl_shape)

    coefficients = ufl.algorithms.extract_coefficients(expression)
    coefficient_numbering = {}
    for i, coeff in enumerate(coefficients):
        coefficient_numbering[coeff] = i

    # Add coefficient numbering to IR
    ir["coefficient_numbering"] = coefficient_numbering

    original_coefficient_positions = []
    original_coefficients = ufl.algorithms.extract_coefficients(original_expression)
    for coeff in coefficients:
        original_coefficient_positions.append(original_coefficients.index(coeff))

    ir["coefficient_names"] = [
        object_names.get(id(obj), f"w{j}") for j, obj in enumerate(coefficients)
    ]

    ir["constant_names"] = [
        object_names.get(id(obj), f"c{j}")
        for j, obj in enumerate(ufl.algorithms.analysis.extract_constants(expression))
    ]

    fs = {}
    for function in tuple(original_coefficients) + tuple(arguments):
        name = object_names.get(id(function), str(function))
        if not str(name).isidentifier():
            raise ValueError(f'Function name "{name}" must be a valid object identifier.')
        el = function.ufl_function_space().ufl_element()
        space = function.ufl_function_space()
        domain = space.ufl_domain()
        cmap = domain.ufl_coordinate_element()
        family = cmap.family_name
        degree = cmap.degree
        value_shape = space.value_shape
        fs[name] = (
            finite_element_hashes[el],
            family,
            degree,
            cmap.cell_type,
            cmap.lagrange_variant,
            value_shape,
        )

    expression_name = object_names.get(id(original_expression), index)

    ir["function_spaces"] = fs
    ir["name_from_uflfile"] = f"expression_{prefix}_{expression_name}"

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
    ir["coefficient_offsets"] = offsets

    ir["integral_type"] = "expression"
    if cell is not None:
        if (tdim := cell.topological_dimension()) == (pdim := points.shape[1]):
            ir["entitytype"] = "cell"
        elif tdim - 1 == pdim:
            ir["entitytype"] = "facet"
        else:
            raise ValueError(
                f"Expression on domain with topological dimension {tdim}"
                + f"with points of dimension {pdim} not supported."
            )
    else:
        # For spatially invariant expressions, all expressions are evaluated in the cell
        ir["entitytype"] = "cell"

    # Build offsets for Constants
    original_constant_offsets = {}
    _offset = 0
    for constant in ufl.algorithms.analysis.extract_constants(original_expression):
        original_constant_offsets[constant] = _offset
        _offset += np.prod(constant.ufl_shape, dtype=int)

    ir["original_constant_offsets"] = original_constant_offsets

    ir["points"] = points

    weights = np.array([1.0] * points.shape[0])
    rule = QuadratureRule(points, weights)
    integrands = {rule: expression}

    if cell is None:
        assert (
            len(ir["original_coefficient_positions"]) == 0
            and len(ir["original_constant_offsets"]) == 0
        )

    expression_ir = compute_integral_ir(
        cell, ir["integral_type"], ir["entitytype"], integrands, tensor_shape, options, visualise
    )
    ir["options"] = options
    ir.update(expression_ir)

    return ExpressionIR(**ir)
