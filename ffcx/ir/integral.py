# Copyright (C) 2013-2025 Martin Sandve Alnæs, Michal Habera and Jørgen S. Dokken
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Main algorithm for building the integral intermediate representation."""

import collections
import logging
import typing
from enum import Enum

import basix
import numpy as np
import numpy.typing as npt
import ufl
from ufl.algorithms.balancing import balance_modifiers
from ufl.checks import is_cellwise_constant
from ufl.classes import QuadratureWeight

from ffcx.definitions import entity_types, supported_integral_types
from ffcx.ir.analysis.factorization import compute_argument_factorization
from ffcx.ir.analysis.graph import ExpressionGraph, build_scalar_graph
from ffcx.ir.analysis.modified_terminals import (
    ModifiedTerminal,
    analyse_modified_terminal,
    is_modified_terminal,
)
from ffcx.ir.analysis.visualise import visualise_graph
from ffcx.ir.elementtables import (
    UniqueTableReferenceT,
    build_optimized_tables,
)
from ffcx.ir.elementtables import (
    table_types as _table_types,
)
from ffcx.ir.representationutils import QuadratureRule
from ffcx.ir.sumfactorisation import (
    SumFactorisationPlan,
    evaluation_stage_sizes,
    factorises,
    plan_sum_factorisation,
)

logger = logging.getLogger("ffcx")


class TensorPart(Enum):
    """What part of the tensor to assemble."""

    full = 1
    diagonal = 4

    @classmethod
    def from_str(cls, value: str):
        """Workaround for string enum prior to Python 3.11."""
        if value == "full":
            return cls.full
        elif value == "diagonal":
            return cls.diagonal
        else:
            raise ValueError(f"Unknown tensor part: {value}")

    def __str__(self) -> str:
        """Get string value."""
        if self == TensorPart.full:
            return "full"
        elif self == TensorPart.diagonal:
            return "diagonal"
        else:
            raise ValueError(f"Unknown cell type: {self}")

    def __int__(self) -> int:
        """Get integer value."""
        return self.value


class ModifiedArgumentDataT(typing.NamedTuple):
    """Modified argument data."""

    ma_index: int
    tabledata: UniqueTableReferenceT


class BlockDataT(typing.NamedTuple):
    """Block data."""

    ttypes: tuple[str, ...]  # list of table types for each block rank
    factor_indices_comp_indices: list[tuple[int, int]]  # list of (factor index, component index)
    all_factors_piecewise: bool  # True if all factors for this block are piecewise
    unames: tuple[str, ...]  # list of unique FE table names for each block rank
    restrictions: tuple[str | None, ...]  # restriction "+" | "-" | None for each block rank
    transposed: bool  # block is the transpose of another
    is_uniform: bool
    ma_data: tuple[ModifiedArgumentDataT, ...]  # used in "full", "safe" and "partial"
    is_permuted: bool  # Do quad points on facets need to be permuted?


class IntermediateIntegrandIR(typing.TypedDict):
    """Intermediate IR for integrands."""

    factorization: ExpressionGraph
    modified_arguments: list[ModifiedTerminal]
    block_contributions: dict[tuple[tuple[int, ...], ...], list[BlockDataT]]
    # How to contract the blocks one direction at a time, if they are
    # sum factorised
    sum_factorisation: SumFactorisationPlan | None


class IntermediateIntegralIR(typing.TypedDict):
    """Intermediate IR for integrals."""

    needs_facet_permutations: bool
    unique_tables: dict[basix.CellType, dict[str, npt.NDArray[np.float64]]]
    unique_table_types: dict[basix.CellType, dict[str, _table_types]]
    integrand: dict[tuple[basix.CellType, QuadratureRule], IntermediateIntegrandIR]


class CommonExpressionIR(typing.NamedTuple):
    """Common-ground for IntegralIR and ExpressionIR."""

    integral_type: supported_integral_types | typing.Literal["expression"]
    entity_type: entity_types  # The entity IR is over
    tensor_shape: list[int]  # The shape of the locally assembled tensor

    # When an expression is processed in UFL, coefficients might drop out of the
    # variational form due to for instance differentiation.
    # This maps from the original coefficient to its new index in the form.
    coefficient_numbering: dict[ufl.Coefficient, int]
    # Each coefficient is accessed from a given offset in the input array
    coefficient_offsets: dict[ufl.Coefficient, int]
    original_constant_offsets: dict[ufl.Constant, int]
    unique_tables: dict[basix.CellType, npt.NDArray[np.float64]]
    unique_table_types: dict[basix.CellType, dict[str, str]]
    integrand: dict[tuple[basix.CellType, QuadratureRule], IntermediateIntegrandIR]
    name: str
    needs_facet_permutations: bool
    shape: list[int]
    coordinate_element_hash: str
    number_coordinate_dofs: int
    # The integration domain's own coordinate element, used to
    # get the coordinate closure dofs for each mixed-dimensional
    # submesh. Avoids having to pack duplicate coordinates.
    integration_domain_coordinate_element: basix.ufl._ElementBase | None = None


def _compute_integral_ir(
    expression: ufl.core.expr.Expr,
    existing_tables: dict[str, npt.NDArray[np.float64]],
    quadrature_rule: QuadratureRule,
    cell: ufl.Cell,
    integral_type: supported_integral_types,
    entity_type: entity_types,
    argument_shape: tuple[int, ...],
    visualise: str | None,
    graph_name: str,
    p: dict,
) -> tuple[
    dict[str, npt.NDArray[np.float64]],
    dict[str, _table_types],
    ExpressionGraph,
    list[int],
    dict[tuple[tuple[int, ...], ...], list[BlockDataT]],
    SumFactorisationPlan | None,
    dict[int, ModifiedTerminal],
    bool,
]:
    # Remove QuadratureWeight terminals from expression and replace with 1.0
    expression = replace_quadratureweight(expression)

    # Build initial scalar list-based graph representation
    S = build_scalar_graph(expression)

    # Build terminal_data from V here before factorization. Then we
    # can use it to derive table properties for all modified
    # terminals, and then use that to rebuild the scalar graph more
    # efficiently before argument factorization. We can build
    # terminal_data again after factorization if that's necessary.

    initial_terminals: dict[int, ModifiedTerminal] = {
        i: analyse_modified_terminal(v["expression"])
        for i, v in S.nodes.items()
        if is_modified_terminal(v["expression"])
    }

    # Check if we have a mixed-dimensional integral.
    # We allow for constants to be defined on other meshes, as they require no mapping.
    is_mixed_dim = False
    coeffs_and_arguments = ufl.algorithms.analysis.extract_arguments(
        expression
    ) + ufl.algorithms.analysis.extract_coefficients(expression)
    domains = map(ufl.domain.extract_unique_domain, coeffs_and_arguments)
    for domain in domains:
        assert domain is not None
        if domain.topological_dimension != cell.topological_dimension:
            is_mixed_dim = True
    # Check geometric terminals for mixed-dimensionality
    for mt in initial_terminals.values():
        if type(mt.terminal) not in (ufl.geometry.SpatialCoordinate, ufl.geometry.Jacobian):
            continue
        domain = ufl.domain.extract_unique_domain(mt.terminal)
        assert domain is not None
        if domain.topological_dimension != cell.topological_dimension:
            is_mixed_dim = True

    mt_table_reference = build_optimized_tables(
        quadrature_rule,
        cell,
        integral_type,
        entity_type,
        initial_terminals.values(),
        existing_tables,
        use_sum_factorization=p["sum_factorization"],
        is_mixed_dim=is_mixed_dim,
        rtol=p["table_rtol"],
        atol=p["table_atol"],
    )
    # Fetch unique tables for this quadrature rule
    table_types: dict[str, _table_types] = {v.name: v.ttype for v in mt_table_reference.values()}
    tables: dict[str, npt.NDArray[np.float64]] = {
        v.name: v.values for v in mt_table_reference.values()
    }

    S_targets = [i for i, v in S.nodes.items() if v.get("target", False)]
    num_components = np.int32(np.prod(expression.ufl_shape))

    if "zeros" in table_types.values():
        # If there are any 'zero' tables, replace symbolically and rebuild graph
        for i, mt in initial_terminals.items():
            # Set modified terminals with zero tables to zero
            tr = mt_table_reference.get(mt)
            if tr is not None and tr.ttype == "zeros":
                S.nodes[i]["expression"] = ufl.as_ufl(0.0)

        # Propagate expression changes using dependency list
        for i, v in S.nodes.items():
            deps = [S.nodes[j]["expression"] for j in S.out_edges[i]]
            if deps:
                v["expression"] = v["expression"]._ufl_expr_reconstruct_(*deps)

        # Recreate expression with correct ufl_shape
        expressions = np.empty(num_components, dtype=CommonExpressionIR)

        for target in S_targets:
            for comp in S.nodes[target]["component"]:
                assert expressions[comp] is None
                expressions[comp] = S.nodes[target]["expression"]
        expression = ufl.as_tensor(expressions.reshape(expression.ufl_shape))
        assert all([expr is not None for expr in expressions])

        # Rebuild scalar list-based graph representation
        S = build_scalar_graph(expression)

    # Output diagnostic graph as pdf
    if visualise:
        visualise_graph(S, f"{graph_name}_S.{visualise}")

    # Compute factorization of arguments
    if p["part"] == TensorPart.diagonal:
        assert len(argument_shape) == 2, "Can only diagonalize bi-linear forms."
        rank = 1
    else:
        rank = len(argument_shape)
    F = compute_argument_factorization(S, rank)

    # Get the 'target' nodes that are factors of arguments, and insert in dict
    FV_targets = [i for i, v in F.nodes.items() if v.get("target", False)]
    argument_factorization: dict[tuple[int, ...], list[tuple[int, int]]] = {}

    for fi in FV_targets:
        # Number of blocks using this factor must agree with number of components
        # to which this factor contributes. I.e. there are more blocks iff there are more
        # components
        assert len(F.nodes[fi]["target"]) == len(F.nodes[fi]["component"])

        k = 0
        for w in F.nodes[fi]["target"]:
            comp = F.nodes[fi]["component"][k]
            argument_factorization[w] = argument_factorization.get(w, [])

            # Store tuple of (factor index, component index)
            argument_factorization[w].append((fi, comp))
            k += 1

    # Get list of indices in F which are the arguments (should be at start)
    _argkeys: set[int] = set()
    for w in argument_factorization:
        _argkeys = _argkeys | set(w)
    argkeys = list(_argkeys)

    # Build set of modified_terminals for each mt factorized vertex in F
    # and attach tables, if appropriate
    for i, v in F.nodes.items():
        expr = v["expression"]
        if is_modified_terminal(expr):
            mt = analyse_modified_terminal(expr)
            F.nodes[i]["mt"] = mt
            tr = mt_table_reference.get(mt)
            if tr is not None:
                F.nodes[i]["tr"] = tr

    # Sum factorise if the rule is a tensor product rule and every block's
    # arguments factorise over its directions. Only the tables of the
    # arguments, and of the terminals then evaluated on the whole rule, keep
    # their tensor factors; other tables are accessed whole.
    dim = len(quadrature_rule.tensor_factors) if quadrature_rule.has_tensor_factors else 0
    sum_factorise = (
        dim > 0
        and len(argument_factorization) > 0
        and all(
            len(ma_indices) > 0
            and len(fi_ci) == 1
            and all(factorises(F.nodes[ai]["tr"], dim) for ai in ma_indices)
            for ma_indices, fi_ci in argument_factorization.items()
        )
    )
    _strip_tensor_factors(
        F,
        lambda i, tr: (
            sum_factorise and (i in _argkeys or (factorises(tr, dim) and not tr.is_piecewise))
        ),
    )

    # Attach 'status' to each node: 'inactive', 'piecewise' or 'varying'
    analyse_dependencies(F, mt_table_reference)

    # Output diagnostic graph as pdf
    if visualise:
        visualise_graph(F, f"{graph_name}_F.{visualise}")

    block_contributions = _block_contributions(F, argument_factorization, tables, rank, p)

    plan = None
    if sum_factorise:
        max_scalars = _max_sum_factorisation_scalars(F, _argkeys, quadrature_rule, p)
        plan = plan_sum_factorisation(quadrature_rule, block_contributions, max_scalars)
        if plan is None:
            # The arrays of the contraction don't fit the stack limit
            _strip_tensor_factors(F)
            block_contributions = _block_contributions(F, argument_factorization, tables, rank, p)

    # Figure out which table names are referenced
    active_table_names = set()
    for i, v in F.nodes.items():
        tr = v.get("tr")
        if tr is not None and F.nodes[i]["status"] != "inactive":
            if tr.has_tensor_factorisation:
                assert tr.tensor_factors is not None
                for t in tr.tensor_factors:
                    active_table_names.add(t.name)
            else:
                active_table_names.add(tr.name)

    # Figure out which table names are referenced in blocks. A sum-factorised
    # contraction reads weighted tables instead of most argument tables.
    if plan is not None:
        active_table_names |= plan.factors
    else:
        for contributions in block_contributions.values():
            for blockdata in contributions:
                active_table_names.update(mad.tabledata.name for mad in blockdata.ma_data)

    active_tables: dict[str, npt.NDArray[np.float64]] = {}
    active_table_types: dict[str, _table_types] = {}

    for name in active_table_names:
        # Drop tables not referenced from modified terminals
        if table_types[name] not in ("zeros", "ones"):
            active_tables[name] = tables[name]
            active_table_types[name] = table_types[name]
    return (
        active_tables,
        active_table_types,
        F,
        argkeys,
        block_contributions,
        plan,
        initial_terminals,
        is_mixed_dim,
    )


def _strip_tensor_factors(F: ExpressionGraph, keep=lambda i, tr: False):
    """Drop the tensor factors of the tables of nodes that are accessed whole."""
    for i, v in F.nodes.items():
        tr = v.get("tr")
        if tr is not None and tr.has_tensor_factorisation and not keep(i, tr):
            v["tr"] = tr._replace(tensor_factors=None)


def _block_contributions(
    F: ExpressionGraph,
    argument_factorization: dict[tuple[int, ...], list[tuple[int, int]]],
    tables: dict[str, npt.NDArray[np.float64]],
    rank: int,
    p: dict,
) -> dict[tuple[tuple[int, ...], ...], list[BlockDataT]]:
    """Return the blocks of the element tensor each factorisation term contributes to."""
    block_contributions: dict[tuple[tuple[int, ...], ...], list[BlockDataT]] = (
        collections.defaultdict(list)
    )
    for ma_indices, fi_ci in sorted(argument_factorization.items()):
        # Get a bunch of information about this term
        if TensorPart.from_str(p["part"]) != TensorPart.diagonal:
            assert rank == len(ma_indices)

        trs = tuple(F.nodes[ai]["tr"] for ai in ma_indices)
        unames = tuple(tr.name for tr in trs)
        ttypes = tuple(tr.ttype for tr in trs)
        assert not any(tt == "zeros" for tt in ttypes)

        _blockmap: list[tuple[int, ...]] = []
        for tr in trs:
            assert tr is not None
            begin = tr.offset
            assert begin is not None
            num_dofs = tr.values.shape[3]
            assert tr.block_size is not None
            dofmap = tuple(begin + i * tr.block_size for i in range(num_dofs))
            _blockmap.append(dofmap)
        blockmap = tuple(_blockmap)

        block_is_uniform = all(tr.is_uniform for tr in trs)

        # Collect relevant restrictions to identify blocks correctly
        # in interior facet integrals
        _block_restrictions: list[str] = []
        for i, ai in enumerate(ma_indices):
            if not trs[i].is_uniform:
                r = F.nodes[ai]["mt"].restriction
                _block_restrictions.append(r)

        block_restrictions: tuple[str, ...] = tuple(_block_restrictions)

        # Check if each *each* factor corresponding to this argument is piecewise
        all_factors_piecewise = all(F.nodes[ifi[0]]["status"] == "piecewise" for ifi in fi_ci)
        block_is_permuted = False
        for name in unames:
            if tables[name].shape[0] > 1:
                block_is_permuted = True
        ma_data = []
        for i, ma in enumerate(ma_indices):
            ma_data.append(ModifiedArgumentDataT(ma, trs[i]))

        block_is_transposed = False  # FIXME: Handle transposes for these block types
        block_unames = unames
        blockdata = BlockDataT(
            ttypes,
            fi_ci,
            all_factors_piecewise,
            block_unames,
            block_restrictions,
            block_is_transposed,
            block_is_uniform,
            tuple(ma_data),
            block_is_permuted,
        )

        # Insert in expr_ir for this quadrature loop
        block_contributions[blockmap].append(blockdata)

    return block_contributions


def _max_sum_factorisation_scalars(
    F: ExpressionGraph, argkeys: set[int], quadrature_rule: QuadratureRule, p: dict
) -> int | None:
    """Return how many scalars the arrays of a sum-factorised contraction may take.

    That is the stack limit less the arrays evaluating coefficients and
    geometry on the whole rule, counting stages shared by terminals reading
    the same dofs with the same 1D tables once, or `None` for no limit.
    The stages are counted last direction first at every point, an upper
    bound on those of `FFCXBackendDefinitions._sum_factorised_lincomb`,
    which contracts point-independent directions first and stores them
    at one point.
    """
    limit = int(p["sum_factorization_stack_limit"])
    if limit <= 0:
        return None
    scalar = np.dtype(p["scalar_type"])
    real = np.finfo(scalar).dtype
    num_points = quadrature_rule.num_points
    evaluation = 0
    stages: set[tuple] = set()
    for i, v in F.nodes.items():
        tr = v.get("tr")
        if i in argkeys or tr is None or not tr.has_tensor_factorisation:
            continue
        if v["status"] != "varying":
            continue
        assert tr.tensor_factors is not None
        terminal = v["mt"].terminal
        num_dofs = [int(factor.values.shape[-1]) for factor in tr.tensor_factors]
        sizes = evaluation_stage_sizes(num_dofs, num_points)
        if isinstance(terminal, ufl.Coefficient):
            itemsize = scalar.itemsize
            key: tuple = (terminal, tr.offset, tr.block_size)
        else:
            # The spatial coordinate and the Jacobian read the coordinate dofs
            itemsize = real.itemsize
            key = (ufl.domain.extract_unique_domain(terminal), tr.offset)
        for d in reversed(range(len(sizes))):
            key += (tr.tensor_factors[d].name,)
            if key not in stages:
                stages.add(key)
                evaluation += sizes[d] * itemsize
    return max(limit - evaluation, 0) // int(scalar.itemsize)


def _graph_name(name, integral_domain, quadrature_rule) -> str:
    """Build the file name stem for an integrand's graph visualisation.

    An integral may be split over several cell types and quadrature
    rules, each with its own graph, so the stem identifies which one it
    is. The quadrature rule identifier is the one used to name tables in
    the generated code, so the two can be matched up.
    """
    parts = [name]
    if integral_domain != "":
        parts.append(integral_domain.name)
    parts.append(f"Q{quadrature_rule.id()}")
    return "_".join(parts)


def compute_integral_ir(
    cell: ufl.Cell,
    integral_type: supported_integral_types,
    entity_type: entity_types,
    integrands: dict[basix.CellType, dict[QuadratureRule, ufl.core.expr.Expr]],
    argument_shape: tuple[int, ...],
    p: dict,
    visualise: str | None,
    name: str,
) -> IntermediateIntegralIR:
    """Compute intermediate representation for an integral.

    Args:
        cell: Cell of integration domain
        integral_type: Type of integral over cell
        entity_type: Corresponding entity of the cell that the integral is over
        integrands: Dictionary mapping a quadrature rule to a sequence of integrands
        argument_shape: Shape of the output tensor of the integral (used for tensor factorization)
        p: Parameters used for clamping tables and for activating sum factorization
        visualise: If set, the image format in which to store the graph
            representation of each integrand, e.g. "pdf" or "png"
        name: Name of the integral or expression, used as the stem of the
            graph file names
    """
    # Data that will populate the intermediate representation.
    needs_facet_permutations: bool = False
    unique_tables: dict[basix.CellType, dict[str, npt.NDArray[np.float64]]] = {}
    unique_table_types: dict[basix.CellType, dict[str, _table_types]] = {}
    integrand_map: dict[tuple[basix.CellType, QuadratureRule], IntermediateIntegrandIR] = {}
    for integral_domain, integrands_on_domain in integrands.items():
        unique_tables[integral_domain] = {}
        unique_table_types[integral_domain] = {}
        for quadrature_rule, integrand in integrands_on_domain.items():
            expression = integrand

            # Rebalance order of nested terminal modifiers
            expression = balance_modifiers(expression)

            (
                active_tables,
                active_table_types,
                F,
                argkeys,
                block_contributions,
                plan,
                initial_terminals,
                is_mixed_dim,
            ) = _compute_integral_ir(
                expression,
                unique_tables[integral_domain],
                quadrature_rule,
                cell,
                integral_type,
                entity_type,
                argument_shape,
                visualise,
                _graph_name(name, integral_domain, quadrature_rule),
                p,
            )

            # Add tables and types for this quadrature rule to global tables dict
            unique_tables[integral_domain].update(active_tables)
            unique_table_types[integral_domain].update(active_table_types)

            # Without sum factorisation, use the tensor product rule as an
            # ordinary rule: its points and weights are already flattened
            if plan is None and quadrature_rule.has_tensor_factors:
                quadrature_rule = QuadratureRule(quadrature_rule.points, quadrature_rule.weights)

            # Build IR dict for the given expressions
            # Store final ir for this num_points
            integrand_map[(integral_domain, quadrature_rule)] = IntermediateIntegrandIR(
                factorization=F,
                modified_arguments=[F.nodes[i]["mt"] for i in argkeys],
                block_contributions=block_contributions,
                sum_factorisation=plan,
            )

            restrictions = [i.restriction for i in initial_terminals.values()]
            if not needs_facet_permutations:
                needs_facet_permutations = (
                    "+" in restrictions and "-" in restrictions
                ) or is_mixed_dim

    return IntermediateIntegralIR(
        needs_facet_permutations=needs_facet_permutations,
        unique_tables=unique_tables,
        unique_table_types=unique_table_types,
        integrand=integrand_map,
    )


def analyse_dependencies(F, mt_unique_table_reference):
    """Analyse dependencies.

    Sets 'status' of all nodes to either: 'inactive', 'piecewise' or 'varying'
    Children of 'target' nodes are either 'piecewise' or 'varying'.
    All other nodes are 'inactive'.
    Varying nodes are identified by their tables ('tr'). All their parent
    nodes are also set to 'varying' - any remaining active nodes are 'piecewise'.
    """
    # Set targets, and dependencies to 'active'
    targets = [i for i, v in F.nodes.items() if v.get("target")]
    for _, v in F.nodes.items():
        v["status"] = "inactive"

    while targets:
        s = targets.pop()
        F.nodes[s]["status"] = "active"
        for j in F.out_edges[s]:
            if F.nodes[j]["status"] == "inactive":
                targets.append(j)

    # Build piecewise/varying markers for factorized_vertices
    varying_ttypes = ("varying", "quadrature", "uniform")
    varying_indices = []
    for i, v in F.nodes.items():
        if v.get("mt") is None:
            continue
        tr = v.get("tr")
        if tr is not None:
            ttype = tr.ttype
            # Check if table computations have revealed values varying over points
            if ttype in varying_ttypes:
                varying_indices.append(i)
            else:
                if ttype not in ("fixed", "piecewise", "ones", "zeros"):
                    raise RuntimeError(f"Invalid ttype {ttype}.")

        elif not is_cellwise_constant(v["expression"]):
            raise RuntimeError("Error " + str(tr))
            # Keeping this check to be on the safe side,
            # not sure which cases this will cover (if any)
            # varying_indices.append(i)

    # Set all parents of active varying nodes to 'varying'
    while varying_indices:
        s = varying_indices.pop()
        if F.nodes[s]["status"] == "active":
            F.nodes[s]["status"] = "varying"
            for j in F.in_edges[s]:
                varying_indices.append(j)

    # Any remaining active nodes must be 'piecewise'
    for _, v in F.nodes.items():
        if v["status"] == "active":
            v["status"] = "piecewise"


def replace_quadratureweight(expression):
    """Remove any QuadratureWeight terminals and replace with 1.0."""
    r = []
    for node in ufl.corealg.traversal.unique_pre_traversal(expression):
        if is_modified_terminal(node) and isinstance(node, QuadratureWeight):
            r.append(node)

    replace_map = {q: 1.0 for q in r}
    return ufl.algorithms.replace(expression, replace_map)
