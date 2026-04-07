# Copyright (C) 2009-2026 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
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

if typing.TYPE_CHECKING:
    from ufl.algorithms.formdata import FormData
    from ufl.core.expr import Expr
import warnings

import basix
import numpy as np
import numpy.typing as npt
import ufl.algorithms
from ufl.classes import Integral
from ufl.sorting import sorted_expr_sum

from ffcx import naming
from ffcx.analysis import ProxyCoefficient, UFLData
from ffcx.definitions import entity_types, supported_integral_types
from ffcx.ir.integral import CommonExpressionIR, TensorPart, compute_integral_ir
from ffcx.ir.representationutils import (
    QuadratureRule,
    create_quadrature_points_and_weights,
)

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
    sub_expressions: list[tuple[ProxyCoefficient, str]]
    proxy_coefficient_sizes: list[int]
    proxy_pack_shape: list[tuple[int, ...]]
    coefficients_in_proxy: list[ufl.Coefficient]
    proxy_coefficient_offsets: list[int]


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


def _group_integrands_by_quadrature_rule(
    integrals: list[ufl.Integral],
    argument_elements: tuple[basix.ufl._ElementBase, ...],
    integral_type: supported_integral_types,
    ufl_cell: ufl.Cell,
    sum_factorization: bool,
) -> dict[basix.CellType, dict[QuadratureRule, list[Expr]]]:
    """Group integrands with the same quadrature rule.

    Given a sequence of integrals, group them by the quadrature
    (after explicit computation of it).
    The grouping is nested, meaning that it is first grouped
    by the cell used for the integration (in most cases the integration domain).
    However, for integrals with vertex-quadrature, the celltype is changed
    to be of the correct sub-entity, so that one fetches the correct quadrature points.

    Args:
        integrals: List of itegrals to group
        argument_elements: The finite elements for the arguments in the form.
            Used to get the correct polyset type in the case of tensor-factorization
            of the element ad thus the quadrature rule.
        integral_type: Type of (FFCx) integral that is performed (not UFL integral type)
        ufl_cell: UFL cell for the integration domain.
        sum_factorization: If True use sum factorization.

    Returns:
        A nested dictionary over integrands[cell type of quadrature][quadrature_rule].
    """
    #
    grouped_integrands: dict[basix.CellType, dict[QuadratureRule, list[Expr]]] = {}
    # NOTE: this variable changes throughout the loop
    cell_type = basix_cell_from_string(ufl_cell.cellname)
    use_sum_factorization = sum_factorization and integral_type == "cell"
    for integral in integrals:
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
            weights = np.full(points.shape[0], cell_volume / points.shape[0], dtype=points.dtype)
            rules[cell_type] = (points, weights, None)
        else:
            degree = md["quadrature_degree"]
            points, weights, tensor_factors = create_quadrature_points_and_weights(
                integral_type,
                ufl_cell,
                degree,
                scheme,
                argument_elements,
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
            if rule not in grouped_integrands[cell_type]:
                grouped_integrands[cell_type][rule] = []
            grouped_integrands[cell_type][rule].append(integral.integrand())
    return grouped_integrands


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
        form_names[fd_index] = naming.form_name(fd.original_form, fd_index, prefix)
        for itg_index, itg_data in enumerate(fd.integral_data):
            integral_names[(fd_index, itg_index)] = naming.integral_name(
                fd.original_form,
                itg_data.integral_type,
                fd_index,
                itg_data.subdomain_id,
                prefix,
            )

    expression_names = {
        org_expr: naming.expression_name((org_expr, points), prefix)
        for (expr, points, org_expr) in analysis.expressions
    }
    irs = [
        _compute_integral_ir(
            fd,
            i,
            analysis.element_numbers.keys(),
            integral_names,
            expression_names,
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
    form_data: FormData,
    form_index: int,
    unique_elements: typing.Iterable[basix.ufl._ElementBase],
    integral_names: dict[tuple[int, int], str],
    expression_names: dict[ufl.core.expr.Expr, str],
    options,
    visualise,
) -> list[IntegralIR]:
    """Compute intermediate representation for form integrals.

    Args:
        form_data: Data from UFL analysis of the form
        form_index: Index of form in the sequence of forms.
        unique_elements: Set of unique elements in the form.
        integral_names: Map from `(form_index, integral_index)` to the name of the integral.
        expression_names: Map from original expression to the name of the expression.
            Used for sub-expressions that are coefficients in the integral (internal proxies).
        options: Options for the intermediate representation. 'part': If the full tensor or
            the diagonal of the tensor should be generated. Only valid for bi-linear forms.
            'sum_factorization': If sum factorization should be used. Only has an effect on cell
            integrals. 'table_atol',`table_rtol': Absolute and relative tolerance for clamping table
            values at -1, 0 or 1.
        visualise: If True, store the graph representation of the integrand in a pdf file
            `S.pdf` and `F.pdf`

    Returns:
         A list of intermediate representations for each integral of the Form.
    """
    _entity_types: dict[supported_integral_types, entity_types] = {
        "cell": "cell",
        "exterior_facet": "facet",
        "interior_facet": "facet",
        "vertex": "vertex",
        "ridge": "ridge",
    }

    # Iterate over groups of integrals
    irs = []
    for itg_data_index, itg_data in enumerate(form_data.integral_data):
        logger.info(f"Computing IR for integral in integral group {itg_data_index}")

        # Cast integral type as long as we are using strings and not StringEnums
        integral_type = typing.cast(supported_integral_types, itg_data.integral_type)

        # Determine what kind of integral we are considering, as well checking that each integrand
        # in the integral is associated with a domain of the same topological dimensio.
        entity_type = _entity_types[integral_type]
        ufl_cell = itg_data.domain.ufl_cell()
        cell_type = basix_cell_from_string(ufl_cell.cellname)
        tdim = ufl_cell.topological_dimension
        assert all(tdim == itg.ufl_domain().topological_dimension for itg in itg_data.integrals)

        # Initial population of what will become a CommonExpressionIR
        expression_ir = {
            "integral_type": integral_type,
            "entity_type": entity_type,
            "shape": (),
            "coordinate_element_hash": itg_data.domain.ufl_coordinate_element().basix_hash(),
            "number_coordinate_dofs": itg_data.domain.ufl_coordinate_element().dim,
        }
        # Initial population of what will become the IntegralIR
        ir = {
            "rank": form_data.rank,
            "enabled_coefficients": itg_data.enabled_coefficients,
            "part": TensorPart.from_str(options["part"]),
        }

        # Determine if the form compiler has been asked to diagonalize a
        # bilinear form (by only assembling the diagonal entries, modify rank if True)
        diagonalise = False
        if form_data.rank == 2 and ir["part"] == TensorPart.diagonal:
            diagonalise = True
            ir["rank"] = 1
            assert form_data.argument_elements[0] == form_data.argument_elements[1], (
                "Can only diagonalise forms with identical arguments."
            )

        # Pre-compute the dimension number of dofs in each local element.
        # Used to compute the shape of the locally assembled tensor, as well the
        # coefficient offset
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

        if diagonalise:  # Modify output tensor shape if diagonalizing
            expression_ir["tensor_shape"] = expression_ir["tensor_shape"][:1]

        # Group integrals by quadrature rule by splitting them into separate integrands
        grouped_integrands = _group_integrands_by_quadrature_rule(
            itg_data.integrals,
            form_data.argument_elements,
            integral_type,
            ufl_cell,
            options["sum_factorization"],
        )

        # Sum up all integrands for a given quadrature rule
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

        # Build coefficient numbering for UFCx interface here, to avoid
        # renumbering in UFL and application of replace mapping.
        # We extract a map from the original coefficients to their
        # new index in the form, as well as where in the flattened
        # input coefficient data a coefficient should start.
        coefficient_numbering: dict[ufl.Coefficient, int] = {}
        coefficient_offsets: dict[ufl.Coefficient, int] = {}
        proxy_coefficient_numbering: dict[ProxyCoefficient, int] = {}
        _proxy_coefficient_offsets: dict[ProxyCoefficient, int] = {}
        _offset_c = 0
        _offset_p = 0
        width = 2 if integral_type in ("interior_facet") else 1
        first_proxy = True
        start_proxy = len(form_data.reduced_coefficients)
        for i, (coeff, el) in enumerate(
            zip(form_data.reduced_coefficients, form_data.coefficient_elements)
        ):
            if isinstance(coeff, ProxyCoefficient):
                proxy_coefficient_numbering[coeff] = i
                _proxy_coefficient_offsets[coeff] = _offset_p
                _offset_p += width * element_dimensions[el]
                start_proxy = i if first_proxy else start_proxy
                first_proxy = False
            else:
                coefficient_numbering[coeff] = i
                coefficient_offsets[coeff] = _offset_c
                _offset_c += width * element_dimensions[el]

        # Add Coefficient numbering and offsets to IR
        expression_ir["coefficient_numbering"] = coefficient_numbering
        expression_ir["coefficient_offsets"] = coefficient_offsets
        expression_ir["proxy_coefficient_numbering"] = proxy_coefficient_numbering
        expression_ir["proxy_coefficient_offsets"] = _proxy_coefficient_offsets

        # Similar procedure with constants, but they are not removed
        # as they usually contain very little data.
        original_constant_offsets = {}
        _offset = 0
        for constant in form_data.original_form.constants():
            original_constant_offsets[constant] = _offset
            _offset += np.prod(constant.ufl_shape, dtype=int)

        expression_ir["original_constant_offsets"] = original_constant_offsets

        # Create map from number of quadrature points -> integrand
        integrand_map: dict[basix.CellType | str, dict[QuadratureRule, Expr]] = {
            cell_type: {rule: integral.integrand() for rule, integral in cell_integrals.items()}
            for cell_type, cell_integrals in sorted_integrals.items()
        }

        # Build more specific intermediate representation
        integral_ir = compute_integral_ir(
            itg_data.domain.ufl_cell(),
            integral_type,
            entity_type,
            integrand_map,
            expression_ir["tensor_shape"],
            options,
            visualise,
        )

        # Check if the coefficients of the integral is a proxy coefficient, and if it is enabled.
        proxy_coefficients: list[ufl.Coefficient] = []
        if itg_data.integral_coefficients is not None and itg_data.enabled_coefficients is not None:
            enabled_itg_coeffcients = [
                c
                for (i, c) in enumerate(itg_data.integral_coefficients)
                if itg_data.enabled_coefficients[i]
            ]
            proxy_coefficients = [
                coeff
                for coeff in form_data.reduced_coefficients[start_proxy:]
                if coeff in enabled_itg_coeffcients
            ]

        # Num coefficients required for the form, including generating proxy coefficients
        proxy_coefficient_offsets = [0]
        ir["coefficients_in_proxy"] = []
        ir["sub_expressions"] = []
        for p_coeff in proxy_coefficients:
            assert isinstance(p_coeff, ProxyCoefficient)
            for coeff in ufl.algorithms.extract_coefficients(p_coeff.operand):
                ir["coefficients_in_proxy"].append(coeff)
            proxy_coefficient_offsets.append(len(ir["coefficients_in_proxy"]))
            ir["sub_expressions"].append((p_coeff, expression_names[p_coeff.operand]))
        ir["proxy_coefficient_offsets"] = proxy_coefficient_offsets

        ir["proxy_pack_shape"] = []
        ir["proxy_coefficient_sizes"] = []
        for p_coeff in proxy_coefficients:
            value_size = int(np.prod(p_coeff.ufl_shape))
            ir["proxy_coefficient_sizes"].append(p_coeff.ufl_element().dim)
            ir["proxy_pack_shape"].append(
                (
                    value_size,
                    p_coeff.ufl_function_space().ufl_element().basix_element.points.shape[0],
                )
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

    ir["num_coefficients"] = len(form_data.original_coefficient_positions)

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
    ufcx_integral_types = (
        "cell",
        "exterior_facet",
        "interior_facet",
        "vertex",
        "ridge",
    )
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
    analysis_expr: tuple[ufl.core.expr.Expr, np.ndarray, ufl.core.expr.Expr],
    index: int,
    prefix: str,
    analysis: UFLData,
    options: dict,
    visualise: bool,
    object_names,
):
    """Compute intermediate representation of an Expression."""
    logger.info(f"Computing IR for Expression {index}")

    # Compute representation
    ir: dict[str, typing.Any] = {}
    base_ir: dict[str, typing.Any] = {}
    original_input_expr = (analysis_expr[2], analysis_expr[1])

    base_ir["name"] = naming.expression_name(original_input_expr, prefix)

    original_expr = analysis_expr[2]
    points = analysis_expr[1]
    expr = analysis_expr[0]

    domain = ufl.domain.extract_unique_domain(expr)
    # This case corresponds to a spatially constant expression
    # without any dependencies
    if domain is None:
        cell = None
    else:
        assert hasattr(domain, "ufl_cell"), "Expression domain does not have a cell."
        cell = domain.ufl_cell()

    # Prepare dimensions of all unique element in expression, including
    # elements for arguments, coefficients and coordinate mappings
    element_dimensions = {
        element: element.dim + element.num_global_support_dofs
        for element in analysis.unique_elements
    }

    # Extract dimensions for elements of arguments only
    arguments = ufl.algorithms.extract_arguments(expr)
    argument_elements = tuple(f.ufl_function_space().ufl_element() for f in arguments)
    argument_dimensions = tuple(element_dimensions[element] for element in argument_elements)

    tensor_shape = argument_dimensions
    base_ir["tensor_shape"] = tensor_shape

    base_ir["shape"] = list(expr.ufl_shape)

    ir["constant_names"] = [
        object_names.get(id(obj), f"c{j}")
        for j, obj in enumerate(ufl.algorithms.analysis.extract_constants(expr))
    ]

    expr_name = object_names.get(id(original_expr), index)
    ir["name_from_uflfile"] = f"expression_{prefix}_{expr_name}"

    if len(argument_elements) > 1:
        raise RuntimeError("Expression with more than one Argument not implemented.")

    coefficients = ufl.algorithms.extract_coefficients(expr)
    original_coefficients = ufl.algorithms.extract_coefficients(original_expr)

    coefficient_numbering: dict[ufl.Coefficient, int] = {}
    coefficient_offsets: dict[ufl.Coefficient, int] = {}
    proxy_coefficient_numbering: dict[ProxyCoefficient, int] = {}
    proxy_coefficient_offsets: dict[ProxyCoefficient, int] = {}
    _offset_c = 0
    _offset_p = 0
    original_coefficient_positions = []
    for i, coeff in enumerate(coefficients):
        el = coeff.ufl_element()
        if isinstance(coeff, ProxyCoefficient):
            proxy_coefficient_numbering[coeff] = i
            proxy_coefficient_offsets[coeff] = _offset_p
            _offset_p += element_dimensions[el]
        else:
            original_coefficient_positions.append(original_coefficients.index(coeff))
            coefficient_numbering[coeff] = i
            coefficient_offsets[coeff] = _offset_c
            _offset_c += element_dimensions[el]

    ir["original_coefficient_positions"] = original_coefficient_positions
    base_ir["coefficient_numbering"] = coefficient_numbering
    base_ir["coefficient_offsets"] = coefficient_offsets
    base_ir["proxy_coefficient_numbering"] = proxy_coefficient_numbering
    base_ir["proxy_coefficient_offsets"] = proxy_coefficient_offsets

    base_ir["integral_type"] = "expression"
    ir["coefficient_names"] = [
        object_names.get(id(obj), f"w{j}") for j, obj in enumerate(coefficients)
    ]

    if cell is not None:
        if (tdim := cell.topological_dimension) == (pdim := points.shape[1]):
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
    assert isinstance(expr_domain, ufl.Mesh) or expr_domain is None, (
        "Expression domain must be a Mesh or None."
    )
    base_ir["coordinate_element_hash"] = (
        expr_domain.ufl_coordinate_element().basix_hash() if expr_domain is not None else 0
    )
    base_ir["number_coordinate_dofs"] = (
        0 if expr_domain is None else expr_domain.ufl_coordinate_element().dim
    )

    weights = np.array([1.0] * points.shape[0])
    rule = QuadratureRule(points, weights)
    integrands: dict[basix.CellType | str, dict[QuadratureRule, ufl.core.expr.Expr]] = {
        "": {rule: expr}
    }

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
