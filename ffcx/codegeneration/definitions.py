# Copyright (C) 2011-2023 Martin Sandve Alnæs, Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCx/UFCx specific variable definitions."""

import itertools
import logging
import math
from collections.abc import Sequence

import numpy as np
import ufl

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration import geometry
from ffcx.codegeneration.optimizer import expr_key, unroll_loop
from ffcx.definitions import entity_types
from ffcx.ir.analysis.modified_terminals import ModifiedTerminal
from ffcx.ir.elementtables import UniqueTableReferenceT
from ffcx.ir.representationutils import QuadratureRule
from ffcx.ir.sumfactorisation import evaluation_stage_sizes

logger = logging.getLogger("ffcx")


# Cap on unrolling, so a curved high-order simplex doesn't blow up compile
# time. Tensor-product geometry always keeps the loop path regardless.
_MAX_UNROLLED_COORDINATE_DOFS = 16
_MAX_UNROLLED_DENSE_COEFFICIENT_DOFS = 8
_MAX_UNROLLED_SPARSE_COEFFICIENT_DOFS = 20
# Expose small 1D reductions and at most one Q2 hex's worth of stage entries.
_MAX_UNROLLED_FACTOR_SIZE = 3
_MAX_UNROLLED_EVALUATION_SIZE = 27


def _should_unroll_coordinate_dofs(num_dofs: int, has_tensor_factorisation: bool) -> bool:
    """Return whether a coordinate-dof linear combination should be unrolled."""
    return not has_tensor_factorisation and num_dofs <= _MAX_UNROLLED_COORDINATE_DOFS


def _coordinate_basis_groups(values: np.ndarray) -> list[tuple[int, list[tuple[int, float]]]]:
    """Group coordinate-basis columns equal up to sign and tabulation round-off.

    Moving the sign into the coordinate combination lets the C compiler hoist
    coordinate differences out of the quadrature loop. Compare all points,
    entities and permutations, with a tolerance at the scale of tabulation
    round-off rather than the looser table-reuse tolerance.
    """
    groups: list[tuple[int, list[tuple[int, float]]]] = []
    roundoff = 16 * np.finfo(values.dtype).eps
    for column in range(values.shape[-1]):
        data = values[..., column]
        scale = np.max(np.abs(data))
        if scale == 0:
            continue
        for representative, members in groups:
            other = values[..., representative]
            tolerance = roundoff * max(scale, np.max(np.abs(other)))
            if np.max(np.abs(data - other)) <= tolerance:
                members.append((column, 1.0))
                break
            if np.max(np.abs(data + other)) <= tolerance:
                members.append((column, -1.0))
                break
        else:
            groups.append((column, [(column, 1.0)]))
    return groups


# Comments on the sections that evaluate terminals, by terminal type
_TERMINAL_DESCRIPTIONS = {
    ufl.Jacobian: "Differentiate the coordinate mapping to obtain Jacobian components.",
    ufl.SpatialCoordinate: "Evaluate physical coordinates from coordinate degrees of freedom.",
    ufl.Coefficient: "Evaluate coefficient values or derivatives from their degrees of freedom.",
}


def _terminal_section(
    mt: ModifiedTerminal,
    statements: list[L.LNode],
    declarations: Sequence[L.Declaration],
    input: list[L.Symbol],
    output: list[L.Symbol],
    annotations: list[L.Annotation] | None = None,
) -> L.Section:
    """Return a section evaluating a terminal, named after its type."""
    terminal = type(mt.terminal)
    return L.Section(
        terminal.__name__,
        statements,
        declarations,
        input,
        output,
        annotations,
        _TERMINAL_DESCRIPTIONS.get(terminal),
    )


def is_point_independent(factor) -> bool:
    """Check whether a tensor factor table has the same values at every point."""
    return bool(np.all(factor.values == factor.values[:, :, :1, :]))


def create_quadrature_index(quadrature_rule, quadrature_index_symbol):
    """Create a multi index for the quadrature loop."""
    ranges = [0]
    name = quadrature_index_symbol.name
    indices = [L.Symbol(name, dtype=L.DataType.INT)]
    if quadrature_rule:
        ranges = [quadrature_rule.weights.size]
        if quadrature_rule.has_tensor_factors:
            dim = len(quadrature_rule.tensor_factors)
            ranges = [factor[1].size for factor in quadrature_rule.tensor_factors]
            indices = [L.Symbol(name + f"{i}", dtype=L.DataType.INT) for i in range(dim)]

    return L.MultiIndex(indices, ranges)


def create_dof_index(tabledata, dof_index_symbol):
    """Create a multi index for the coefficient dofs."""
    name = dof_index_symbol.name
    if tabledata.has_tensor_factorisation:
        dim = len(tabledata.tensor_factors)
        ranges = [factor.values.shape[-1] for factor in tabledata.tensor_factors]
        indices = [L.Symbol(f"{name}{i}", dtype=L.DataType.INT) for i in range(dim)]
    else:
        ranges = [tabledata.values.shape[-1]]
        indices = [L.Symbol(name, dtype=L.DataType.INT)]

    return L.MultiIndex(indices, ranges)


class FFCXBackendDefinitions:
    """FFCx specific code definitions."""

    entity_type: entity_types

    def __init__(self, entity_type: entity_types, integral_type: str, access, options):
        """Initialise."""
        # Store ir and options
        self.integral_type = integral_type
        self.entity_type = entity_type
        self.access = access
        self.options = options

        # Stages of sum-factorised evaluations, by quadrature rule, dofs and
        # the 1D tables contracted with, shared by the terminals of a kernel
        self._stages: dict[tuple, L.Symbol] = {}

        # called, depending on the first argument type.
        self.handler_lookup = {
            ufl.coefficient.Coefficient: self.coefficient,
            ufl.geometry.Jacobian: self._define_coordinate_dofs_lincomb,
            ufl.geometry.SpatialCoordinate: self.spatial_coordinate,
            ufl.constant.Constant: self.pass_through,
            ufl.geometry.CellVertices: self.pass_through,
            ufl.geometry.FacetEdgeVectors: self.pass_through,
            ufl.geometry.CellEdgeVectors: self.pass_through,
            ufl.geometry.CellFacetJacobian: self.pass_through,
            ufl.geometry.CellRidgeJacobian: self.pass_through,
            ufl.geometry.ReferenceCellVolume: self.pass_through,
            ufl.geometry.ReferenceFacetVolume: self.pass_through,
            ufl.geometry.ReferenceCellEdgeVectors: self.pass_through,
            ufl.geometry.ReferenceFacetEdgeVectors: self.pass_through,
            ufl.geometry.ReferenceNormal: self.pass_through,
            ufl.geometry.CellOrientation: self.pass_through,
            ufl.geometry.FacetOrientation: self.pass_through,
        }

    @property
    def symbols(self):
        """Return formatter."""
        return self.access.symbols

    def get(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> list[L.Section]:
        """Return the sections defining a terminal."""
        # Call appropriate handler, depending on the type of terminal
        terminal = mt.terminal
        ttype = type(terminal)

        # Look for parent class of ttype or direct handler
        while ttype not in self.handler_lookup and ttype.__bases__:
            ttype = ttype.__bases__[0]

        # Get the handler from the lookup, or None if not found
        handler = self.handler_lookup.get(ttype)  # type: ignore

        if handler is None:
            raise NotImplementedError(f"No handler for terminal type: {ttype}")

        # Call the handler
        return handler(mt, tabledata, quadrature_rule, access)  # type: ignore

    def coefficient(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> list[L.Section]:
        """Return definition code for coefficients."""
        # For applying tensor product to coefficients, we need to know
        # if the coefficient has a tensor factorisation and if the
        # quadrature rule has a tensor factorisation. If both are true,
        # we can apply the tensor product to the coefficient.

        iq_symbol = self.symbols.quadrature_loop_index
        ic_symbol = self.symbols.coefficient_dof_sum_index

        iq = create_quadrature_index(quadrature_rule, iq_symbol)
        ic = create_dof_index(tabledata, ic_symbol)

        # Get properties of tables
        ttype = tabledata.ttype
        num_dofs = tabledata.values.shape[3]
        bs = tabledata.block_size
        begin = tabledata.offset
        assert bs is not None
        assert begin is not None
        end = begin + bs * (num_dofs - 1) + 1

        if ttype == "zeros":
            logger.debug("Not expecting zero coefficients to get this far.")
            return []

        # For a constant coefficient we reference the dofs directly, so
        # no definition needed
        if ttype == "ones" and end - begin == 1:
            return []

        assert begin < end

        if tabledata.has_tensor_factorisation:
            dofs = self.symbols.coefficient_dof_access(mt.terminal, begin)
            return self._sum_factorised_lincomb(
                mt, tabledata, quadrature_rule, access, dofs, bs, L.DataType.SCALAR
            )

        # Expand small reductions to expose SIMD across quadrature points,
        # including the evaluation/contraction split for varying geometry.
        # Larger reductions are expanded only when some basis derivatives
        # vanish at every point, or the table is constant so every term
        # folds to a literal. Preserve tensor-product factorisations.
        if (
            not tabledata.has_tensor_factorisation
            and num_dofs <= _MAX_UNROLLED_SPARSE_COEFFICIENT_DOFS
        ):
            active = np.flatnonzero(tabledata.values.any(axis=(0, 1, 2)))
            domain = ufl.domain.extract_unique_domain(mt.terminal)
            assert isinstance(domain, ufl.Mesh)
            unroll_dense = (
                num_dofs <= _MAX_UNROLLED_DENSE_COEFFICIENT_DOFS
                and not domain.is_piecewise_linear_simplex_domain()
            )
            constant = tabledata.is_piecewise and tabledata.is_uniform and not tabledata.is_permuted
            if len(active) < num_dofs or unroll_dense or constant:
                terms = []
                tables = []
                for k in active:
                    ic_k = L.MultiIndex([L.LiteralInt(int(k))], [num_dofs])
                    FE, term_tables = self.access.table_access(
                        tabledata, self.entity_type, mt.restriction, iq, ic_k
                    )
                    dof = self.symbols.coefficient_dof_access(mt.terminal, int(k) * bs + begin)
                    terms.append(L.float_product([dof, FE]))
                    tables.extend(term_tables)
                return self._unrolled_coefficient(mt, access, begin, terms, tables)
        # Get access to element table
        FE, tables = self.access.table_access(tabledata, self.entity_type, mt.restriction, iq, ic)
        dof_access: L.ArrayAccess = self.symbols.coefficient_dof_access(
            mt.terminal, (ic.global_index) * bs + begin
        )

        declaration: list[L.Declaration] = [L.VariableDecl(access, 0.0)]
        body = [L.AssignAdd(access, dof_access * FE)]
        code = [L.create_nested_for_loops([ic], body)]

        input = [dof_access.array, *tables]
        output = [access]
        annotations = [L.Annotation.fuse]

        # assert input and output are Symbol objects
        assert all(isinstance(i, L.Symbol) for i in input)
        assert all(isinstance(o, L.Symbol) for o in output)

        return [_terminal_section(mt, code, declaration, input, output, annotations)]

    def _unrolled_coefficient(
        self,
        mt: ModifiedTerminal,
        access: L.Symbol,
        begin: int,
        terms: list[L.LExpr],
        tables: list[L.Symbol],
    ) -> list[L.Section]:
        """Define a coefficient as an unrolled sum of dof terms."""
        return [
            _terminal_section(
                mt,
                [],
                [L.VariableDecl(access, L.Sum(terms))],
                [
                    self.symbols.coefficient_dof_access(mt.terminal, begin).array,
                    *dict.fromkeys(tables),
                ],
                [access],
            )
        ]

    def _define_coordinate_dofs_lincomb(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> list[L.Section]:
        """Define x or J as a linear combination of coordinate dofs with given table data."""
        # Get properties of domain
        domain = ufl.domain.extract_unique_domain(mt.terminal)
        assert isinstance(domain, ufl.Mesh)
        coordinate_element = domain.ufl_coordinate_element()
        num_scalar_dofs = coordinate_element.sub_elements[0].dim

        num_dofs = tabledata.values.shape[3]
        begin = tabledata.offset
        assert begin is not None

        assert num_scalar_dofs == num_dofs

        # Find table name
        ttype = tabledata.ttype

        assert ttype != "zeros"

        # Get access to element table
        ic_symbol = self.symbols.coefficient_dof_sum_index
        iq_symbol = self.symbols.quadrature_loop_index
        iq = create_quadrature_index(quadrature_rule, iq_symbol)

        dof_access = L.Symbol("coordinate_dofs", dtype=L.DataType.REAL)

        parent_element = self.access.integration_domain_coordinate_element

        # coordinate dofs is always 3d
        dim = 3
        offset = 0
        if mt.restriction == "-":
            # `coordinate_dofs` holds the two cells of the *integration
            # domain* back to back. If terminal is on a submesh,
            # we need the offset to get the second parent cell
            restriction_dofs = (
                num_scalar_dofs if parent_element is None else parent_element.sub_elements[0].dim
            )
            offset = restriction_dofs * dim

        # A submesh's coordinate dofs are a subset of the integration
        # domain's, so gather them through a closure-dofs table rather
        # than indexing `coordinate_dofs` directly.
        closure_table = None
        if (
            parent_element is not None
            and (codim := parent_element.cell.topological_dimension - domain.topological_dimension)
            != 0
        ):
            table_kind, expected_codim = geometry.CLOSURE_DOFS_TABLES.get(
                self.entity_type, (None, None)
            )
            if table_kind is None or codim != expected_codim:
                raise NotImplementedError(
                    "Cannot gather a mixed-dimensional submesh's coordinate dofs: coefficient "
                    f"domain has codimension {codim} relative to the integration domain, but "
                    f"the integral's entity type is {self.entity_type!r}."
                )
            parent_cellname = parent_element.cell.cellname
            closure_table = L.Symbol(f"{parent_cellname}_{table_kind}", dtype=L.DataType.INT)
            entity = self.symbols.entity(self.entity_type, mt.restriction)
            # A vertex table has one row, so the runtime permutation
            # would be meaningless and possibly out of range.
            perm = (
                L.LiteralInt(0)
                if domain.topological_dimension == 0
                else self.symbols.entity_permutation(mt.restriction)
            )
            closure_index = closure_table[perm][entity]

        if closure_table is None and tabledata.has_tensor_factorisation:
            return self._sum_factorised_lincomb(
                mt,
                tabledata,
                quadrature_rule,
                access,
                dof_access[begin + offset],
                dim,
                L.DataType.REAL,
            )

        # Map a submesh-local scalar dof index to its coordinate_dofs index.
        if closure_table is None:

            def _dof(local_index):
                return local_index
        else:

            def _dof(local_index):
                return closure_index[local_index]

        code: list[L.LNode]
        if ttype == "ones":
            # Point meshes have DG-0 basis functions, so the table
            # for this element has been dropped (as it is all ones).
            # A coordinate basis is a partition of unity, so an all-ones
            # table can only come from a single-dof element.
            assert num_dofs == 1
            declaration = [L.VariableDecl(access, dof_access[_dof(0) * dim + begin + offset])]
            code = []
            input = [dof_access]
        elif not _should_unroll_coordinate_dofs(num_dofs, tabledata.has_tensor_factorisation):
            # Many coordinate dofs: keep the runtime loop instead of one
            # literal term per dof.
            ic = create_dof_index(tabledata, ic_symbol)
            FE, tables = self.access.table_access(
                tabledata, self.entity_type, mt.restriction, iq, ic
            )
            code = []
            declaration = [L.VariableDecl(access, 0.0)]
            body = [
                L.AssignAdd(access, dof_access[_dof(ic.global_index) * dim + begin + offset] * FE)
            ]
            code = [L.create_nested_for_loops([ic], body)]
            input = [dof_access, *tables]
        else:
            # Few dofs: emit a single literal-indexed sum instead of a
            # runtime loop. This avoids GCC's vectoriser mis-vectorising a
            # tiny fixed-trip-count reduction into an expensive
            # permute-then-horizontal-sum sequence, and lets the table's
            # exact 0.0/+-1.0 low-order values constant-fold away.
            terms = []
            coord_tables: list[L.Symbol] = []
            # Keep four-coordinate reductions intact: they already fit a
            # short SIMD reduction, and pairing can disrupt that lowering.
            constant_table = (
                tabledata.is_piecewise and tabledata.is_uniform and not tabledata.is_permuted
            )
            groups = (
                [(k, [(k, 1.0)]) for k in range(num_dofs)]
                if constant_table or num_dofs <= 4
                else _coordinate_basis_groups(tabledata.values)
            )
            for k, members in groups:
                ic_k = L.MultiIndex([L.LiteralInt(k)], [num_dofs])
                FE_k, tables_k = self.access.table_access(
                    tabledata, self.entity_type, mt.restriction, iq, ic_k
                )
                if L.is_zero_lexpr(FE_k):
                    continue
                for t in tables_k:
                    if t not in coord_tables:
                        coord_tables.append(t)
                coordinate = dof_access[_dof(k) * dim + begin + offset]
                if len(members) > 1:
                    coordinate = L.Sum(
                        [
                            L.float_product(
                                [dof_access[_dof(j) * dim + begin + offset], L.LiteralFloat(sign)]
                            )
                            for j, sign in members
                        ]
                    )
                terms.append(L.float_product([coordinate, FE_k]))
            declaration = [L.VariableDecl(access, L.Sum(terms))]
            code = []
            input = [dof_access, *coord_tables]

        if closure_table is not None:
            input.append(closure_table)

        output = [access]
        annotations = [L.Annotation.fuse]

        # assert input and output are Symbol objects
        assert all(isinstance(i, L.Symbol) for i in input)
        assert all(isinstance(o, L.Symbol) for o in output)

        return [_terminal_section(mt, code, declaration, input, output, annotations)]

    def _sum_factorised_lincomb(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
        dofs: L.ArrayAccess,
        stride: int,
        dtype: L.DataType,
    ) -> list[L.Section]:
        """Evaluate a linear combination of dofs on the whole tensor product rule.

        The IR keeps the tensor factors of a table only if its terminal is
        evaluated this way.

        Evaluated at each point, a coefficient, or the spatial coordinate or
        Jacobian, sums over all its dofs, costing O(Q^d N^d). Contracting
        its dofs with the 1D tables one
        direction at a time, last direction first, costs O(N^d Q + ... +
        Q^d N); in 2D

            t[i0, q1] = sum_i1 B1[q1, i1] w[i0, i1]
            f[q0, q1] = sum_i0 B0[q0, i0] t[i0, q1]

        Returns a section evaluating f, annotated to run before the
        quadrature loop, and a section reading f at each point.
        """
        iq = create_quadrature_index(quadrature_rule, self.symbols.quadrature_loop_index)
        dim = iq.dim
        assert tabledata.tensor_factors is not None
        num_dofs = [factor.values.shape[-1] for factor in tabledata.tensor_factors]
        # A derivative of a low-order coordinate basis is often constant
        # in its own direction. Store and evaluate that direction once.
        num_points = [
            1 if is_point_independent(factor) else n
            for factor, n in zip(tabledata.tensor_factors, iq.sizes)
        ]
        assert all(
            factor.values.shape[2] == n for factor, n in zip(tabledata.tensor_factors, iq.sizes)
        )

        # Contracting direction d leaves the dofs of directions < d and the
        # points of directions >= d
        sizes = evaluation_stage_sizes(num_dofs, num_points)

        def loop(name: str, extent: int):
            """Return a loop index and a function wrapping a body in its loop.

            A loop of one iteration is left out, its index being zero.
            """
            if extent == 1:
                return L.LiteralInt(0), lambda body: body
            index = L.Symbol(f"sfc_{name}", dtype=L.DataType.INT)
            return index, lambda body: [L.ForRange(index, 0, extent, body)]

        source, offset = dofs.array, dofs.indices[0]
        inputs: list[L.Symbol] = []
        name = f"{access.name}_{quadrature_rule.id()}"
        code: list[L.LNode] = []
        # Terminals reading the same dofs with the same 1D tables, such as x
        # and the derivatives of x, share their first stages
        key: tuple = (quadrature_rule.id(), expr_key(dofs), stride)
        code += [
            L.Comment(f"--- {type(mt.terminal).__name__} component: {access.name} ---"),
            L.Comment(
                "Apply the one-dimensional basis tables to the dofs, "
                "storing point-independent directions only once."
            ),
        ]
        derivative_orders = [mt.local_derivatives.count(d) for d in range(dim)]
        if isinstance(mt.terminal, ufl.geometry.Jacobian):
            row, column = mt.component
            derivative_orders[column] += 1
            code.append(L.Comment(f"Compute J[{row},{column}] = d x{row}/d xi{column}."))
            if mt.local_derivatives:
                code.append(
                    L.Comment(
                        f"Also apply reference derivatives in directions {mt.local_derivatives}."
                    )
                )

        contracted: set[int] = set()

        def stage_comments(direction: int, stage: L.Symbol) -> list[L.LNode]:
            """Describe logical tensor axes before and after this dof contraction."""

            def layout(done: set[int]) -> str:
                axes, shape = [], []
                for d in range(dim):
                    if d in done:
                        if num_points[d] == 1:
                            continue
                        axes.append(f"q{d}")
                        shape.append(num_points[d])
                    else:
                        axes.append(f"i{d}")
                        shape.append(num_dofs[d])
                indices = "".join(f"[{axis}]" for axis in axes) or "scalar"
                return f"{indices}, shape {shape}"

            done = contracted | {direction}
            comments: list[L.LNode] = [
                L.Comment(f"--- Direction {direction}: evaluate {stage.name} ---"),
                L.Comment(f"Logical input layout: {layout(contracted)}."),
                L.Comment(f"Logical output layout: {layout(done)}; last index contiguous."),
            ]
            assert tabledata.tensor_factors is not None
            factor = tabledata.tensor_factors[direction]
            if (
                num_points[direction] == 1
                and num_dofs[direction] == 2
                and np.array_equal(factor.values[0, 0, 0], [-1.0, 1.0])
            ):
                comments.append(
                    L.Comment(
                        f"Subtract the i{direction}=0 slice from the i{direction}=1 slice: "
                        "the Q1 reference derivative has weights [-1, +1]."
                    )
                )
            else:
                operation = (
                    f"Differentiate (order {derivative_orders[direction]})"
                    if derivative_orders[direction]
                    else "Interpolate"
                )
                comments.append(
                    L.Comment(
                        f"{operation} in reference direction {direction}, "
                        f"summing over its {num_dofs[direction]} dofs with the 1D basis table."
                    )
                )
            independent = [f"q{d}" for d in sorted(done) if num_points[d] == 1]
            if independent:
                comments.append(
                    L.Comment(
                        f"Independent of {', '.join(independent)}; omit these axes from storage "
                        "and reuse the values at every point in those directions."
                    )
                )
            return comments

        directions = list(reversed(range(dim)))
        if (
            dim in (2, 3)
            and all(n == 2 for n in num_dofs)
            and all(n <= 3 for n in num_points)
            and any(num_points[d] == 1 for d in range(dim - 1))
        ):
            # Form differences of dofs before interpolating the remaining
            # directions when the derivative table is point-independent.
            result_key = key + tuple(f.name for f in reversed(tabledata.tensor_factors))
            if result_key not in self._stages:
                order = sorted(directions, key=lambda d: num_points[d] != 1)
                shape = list(num_dofs)
                partial_key = key
                for d in order:
                    partial_key += ((d, tabledata.tensor_factors[d].name),)
                    previous = list(shape)
                    shape[d] = num_points[d]
                    if partial_key not in self._stages:
                        if source not in inputs:
                            inputs.append(source)
                        stage = L.Symbol(f"{name}_r{d}", dtype=dtype)
                        code += stage_comments(d, stage)
                        code.append(L.ArrayDecl(stage, sizes=math.prod(shape)))
                        for node_index in itertools.product(*(range(n) for n in shape)):
                            terms = []
                            for dof in range(num_dofs[d]):
                                source_index = list(node_index)
                                source_index[d] = dof
                                flat = L.MultiIndex(source_index, previous).global_index
                                table = tabledata.tensor_factors[d].values[0, 0, node_index[d], dof]
                                terms.append(
                                    L.float_product(
                                        [
                                            source[flat * stride + offset],
                                            L.LiteralFloat(float(table)),
                                        ]
                                    )
                                )
                            target = L.MultiIndex(list(node_index), shape).global_index
                            code.append(L.Assign(stage[target], L.Sum(terms)))
                        self._stages[partial_key] = stage
                    source, stride, offset = self._stages[partial_key], 1, 0
                    contracted.add(d)
                self._stages[result_key] = source
            source = self._stages[result_key]
            directions = []
        for d in directions:
            key += (tabledata.tensor_factors[d].name,)
            if key in self._stages:
                source, stride, offset = self._stages[key], 1, 0
                contracted.add(d)
                continue
            if source not in inputs:
                inputs.append(source)
            stage = L.Symbol(f"{name}_c{d}", dtype=dtype)
            self._stages[key] = stage
            code += stage_comments(d, stage)
            contracted.add(d)
            code.append(L.ArrayDecl(stage, sizes=sizes[d]))
            inner = math.prod(num_points[d + 1 :])
            a, a_loop = loop("a", math.prod(num_dofs[:d]))
            q, q_loop = loop("q", num_points[d])
            width = num_points[-1] if d < dim - 1 else 1
            b, b_loop = loop("b", inner // width)

            unroll = num_dofs[d] <= _MAX_UNROLLED_FACTOR_SIZE
            short_rule = unroll and num_points[d] <= _MAX_UNROLLED_FACTOR_SIZE
            expand_stage = short_rule and sizes[d] <= _MAX_UNROLLED_EVALUATION_SIZE
            points = [L.LiteralInt(k) for k in range(num_points[d])] if short_rule else [q]
            point_code: list[L.LNode] = []
            for point_index in points:
                i = L.Symbol("sfc_i", dtype=L.DataType.INT)
                initial, updates = [], []
                for lane in range(width):
                    point = b * width + lane
                    out = stage[(a * num_points[d] + point_index) * inner + point]
                    terms = []
                    indices = (
                        [L.LiteralInt(k) for k in range(num_dofs[d])]
                        if unroll
                        else [L.LiteralInt(0), i]
                    )
                    for index in indices:
                        value = source[
                            ((a * num_dofs[d] + index) * inner + point) * stride + offset
                        ]
                        if short_rule:
                            assert isinstance(point_index, L.LiteralInt)
                            assert isinstance(index, L.LiteralInt)
                            table = L.LiteralFloat(
                                tabledata.tensor_factors[d].values[
                                    0, 0, int(point_index.value), int(index.value)
                                ]
                            )
                            if L.is_zero_lexpr(table):
                                continue
                        else:
                            table, symbol = self.access.tensor_factor_access(
                                tabledata, d, point_index, index
                            )
                            if symbol not in inputs:
                                inputs.append(symbol)
                        terms.append(L.float_product([value, table]))
                    if unroll:
                        initial.append(L.Assign(out, L.Sum(terms) if terms else 0.0))
                    else:
                        initial.append(L.Assign(out, terms[0]))
                        updates.append(L.AssignAdd(out, terms[1]))

                # Short reductions expose literal table entries and independent
                # output points. Larger stages retain their outer loops.
                body = initial if expand_stage else b_loop(initial)
                if updates and num_dofs[d] > 1:
                    body += [L.ForRange(i, 1, num_dofs[d], b_loop(updates))]
                point_code += body
            stage_code = point_code if short_rule else q_loop(point_code)
            if expand_stage:
                # Literal indices let the compiler replace small temporary
                # arrays by registers instead of vectorising tiny loops.
                for index, extent in ((b, inner // width), (a, math.prod(num_dofs[:d]))):
                    if isinstance(index, L.Symbol):
                        stage_code = unroll_loop(stage_code, index, extent)
            else:
                stage_code = a_loop(stage_code)
            code += stage_code
            source, stride, offset = stage, 1, 0

        terminal = type(mt.terminal).__name__
        evaluation = L.Section(
            f"{terminal} on the tensor product rule",
            code,
            [],
            inputs,
            [source],
            [L.Annotation.factorize],
        )
        read = _terminal_section(
            mt,
            [],
            [
                L.VariableDecl(
                    access,
                    source[
                        L.MultiIndex(
                            [
                                iq.local_index(d) if n > 1 else L.LiteralInt(0)
                                for d, n in enumerate(num_points)
                            ],
                            num_points,
                        )
                    ],
                )
            ],
            [source],
            [access],
            [L.Annotation.fuse],
        )
        return [evaluation, read]

    def spatial_coordinate(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> list[L.Section]:
        """Return definition code for the physical spatial coordinates.

        If physical coordinates are given:
          No definition needed.

        If reference coordinates are given:
          x = sum_k xdof_k xphi_k(X)

        If reference facet coordinates are given:
          x = sum_k xdof_k xphi_k(Xf)
        """
        if self.integral_type in ufl.custom_integral_types:
            # FIXME: Jacobian may need adjustment for custom_integral_types
            if mt.local_derivatives:
                raise ValueError("FIXME: Jacobian in custom integrals is not implemented.")
            return []
        else:
            return self._define_coordinate_dofs_lincomb(mt, tabledata, quadrature_rule, access)

    def jacobian(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> list[L.Section]:
        """Return definition code for the Jacobian of x(X)."""
        return self._define_coordinate_dofs_lincomb(mt, tabledata, quadrature_rule, access)

    def pass_through(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
        access: L.Symbol,
    ) -> list[L.Section]:
        """Return definition code for pass through terminals."""
        return []
