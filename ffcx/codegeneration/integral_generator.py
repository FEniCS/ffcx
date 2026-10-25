# Copyright (C) 2015-2024 Martin Sandve Alnæs, Michal Habera, Igor Baratta, Chris Richardson
#
# Modified by Jørgen S. Dokken, 2024, 2026
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Integral generator."""

import collections
import logging
import math
import platform
import sys
from numbers import Integral
from typing import Any

import basix
import ufl

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration import geometry
from ffcx.codegeneration.definitions import create_dof_index, create_quadrature_index
from ffcx.codegeneration.optimizer import (
    expr_key,
    optimize,
    power_of_two_cse_statements,
    prune_dead_scalars,
    reciprocal_cse_statements,
    scalar_declaration_names,
)
from ffcx.ir.elementtables import piecewise_ttypes
from ffcx.ir.integral import BlockDataT, CommonExpressionIR, TensorPart
from ffcx.ir.representation import IntegralIR
from ffcx.ir.representationutils import QuadratureRule
from ffcx.ir.sumfactorisation import ContractionGroup, SumFactorisationPlan

logger = logging.getLogger("ffcx")


_FW_CACHE_BUDGET_BYTES = 32768
_FW_CACHE_BYTES_PER_SCALAR = 16
_UNSUPPORTED_VECTOR_MATH_MIN_CONTRACTION_WORK = 256
# GCC can lower these calls to glibc libmvec on supported x86-64 Linux
# targets. Other hosts commonly keep scalar libm calls, in which case a
# second loop and an fw cache are pure overhead.
_VECTORIZABLE_MATH_FUNCTIONS = frozenset(
    {
        "acos",
        "acosh",
        "asin",
        "asinh",
        "atan",
        "atanh",
        "cos",
        "cosh",
        "erf",
        "exp",
        "ln",
        "sin",
        "sinh",
        "tan",
        "tanh",
    }
)


def _fw_cache_tile_size(n_fw: int, num_points: int) -> int | None:
    """Return a cache-tile length that stays within the fixed stack budget.

    ``None`` means that even one cached value per intermediate would exceed
    the budget, so the quadrature loop should remain fused.
    """
    assert n_fw > 0
    assert num_points > 0

    max_tile = _FW_CACHE_BUDGET_BYTES // (n_fw * _FW_CACHE_BYTES_PER_SCALAR)
    if max_tile == 0:
        return None

    tile = min(max_tile, num_points)
    # Round down to a SIMD-friendly width when it still fits the budget.
    if tile >= 8:
        tile = tile // 8 * 8
    return tile


def extract_dtype(v, vops: list[Any]):
    """Extract dtype from ufl expression v and its operands."""
    dtypes = []
    for op in vops:
        if hasattr(op, "dtype"):
            dtypes.append(op.dtype)
        elif hasattr(op, "symbol"):
            dtypes.append(op.symbol.dtype)
        elif isinstance(op, Integral):
            dtypes.append(L.DataType.INT)
        else:
            raise RuntimeError(f"Not expecting this type of operand {type(op)}")
    is_cond = isinstance(v, ufl.classes.Condition)
    if is_cond:
        return L.DataType.BOOL
    is_real = isinstance(v, (ufl.classes.Real, ufl.classes.Imag))
    if is_real:
        return L.DataType.REAL
    return L.merge_dtypes(dtypes)


def _contains_vectorizable_mathfunction(node) -> bool:
    """Return whether *node* contains a math call targeted by libmvec."""
    if isinstance(node, list):
        return any(_contains_vectorizable_mathfunction(n) for n in node)
    if isinstance(node, L.MathFunction):
        return node.function in _VECTORIZABLE_MATH_FUNCTIONS or any(
            _contains_vectorizable_mathfunction(arg) for arg in node.args
        )
    if isinstance(node, (L.StatementList, L.Section)):
        return _contains_vectorizable_mathfunction(node.statements)
    if isinstance(node, L.ForRange):
        return _contains_vectorizable_mathfunction(node.body)
    if isinstance(node, L.VariableDecl):
        return node.value is not None and _contains_vectorizable_mathfunction(node.value)
    if isinstance(node, (L.ArrayDecl, L.Comment)):
        return False
    if isinstance(node, L.Statement):
        return _contains_vectorizable_mathfunction(node.expr)
    if isinstance(node, L.NaryOp):
        return any(_contains_vectorizable_mathfunction(a) for a in node.args)
    if isinstance(node, L.BinOp):
        return _contains_vectorizable_mathfunction(node.lhs) or _contains_vectorizable_mathfunction(
            node.rhs
        )
    if isinstance(node, L.PrefixUnaryOp):
        return _contains_vectorizable_mathfunction(node.arg)
    return False


def _supports_vector_math_loop_split() -> bool:
    """Return whether the current JIT host has the required vector libm ABI."""
    return sys.platform.startswith("linux") and platform.machine().lower() in {"amd64", "x86_64"}


def _contraction_work(node) -> int:
    """Estimate statically repeated scalar statements in a contraction tree."""
    if isinstance(node, list):
        return sum(_contraction_work(child) for child in node)
    if isinstance(node, (L.Section, L.StatementList)):
        return _contraction_work(node.statements)
    if isinstance(node, L.ForRange):
        if not isinstance(node.begin, L.LiteralInt) or not isinstance(node.end, L.LiteralInt):
            return 0
        return int(max(0, node.end.value - node.begin.value)) * _contraction_work(node.body)
    return int(isinstance(node, L.Statement))


class IntegralGenerator:
    """Integral generator."""

    def __init__(self, ir: IntegralIR, backend):
        """Initialise."""
        # Store ir
        self.ir = ir

        # Backend specific plugin with attributes
        # - symbols: for translating ufl operators to target language
        # - definitions: for defining backend specific variables
        # - access: for accessing backend specific variables
        self.backend = backend

        # Set of operator names code has been generated for, used in the
        # end for selecting necessary includes
        self._ufl_names: set[str] = set()

        # Initialize lookup tables for variable scopes
        self.init_scopes()

        # Cache
        self.temp_symbols: dict[Any, L.Symbol] = {}

        # Piecewise symbol name -> (lhs, rhs) of its defining two-factor
        # product, if any. Lets generate_block_parts spot several block
        # entries sharing one factor (e.g. several metric-tensor entries all
        # scaled by |det J|) and fuse it with the quadrature weight once.
        self._piecewise_mul_operands: dict[str, tuple[L.LExpr, L.LExpr]] = {}

        # Set of counters used for assigning names to intermediate
        # variables
        self.symbol_counters: dict[str, int] = collections.defaultdict(int)

    def init_scopes(self):
        """Initialize variable scope dicts."""
        # Reset variables, separate sets for each quadrature rule
        self.scopes = {
            quadrature_rule: {} for quadrature_rule in self.ir.expression.integrand.keys()
        }
        self.scopes[(None, None)] = {}

    def set_var(self, quadrature_rule, domain, v, vaccess):
        """Set a new variable in variable scope dicts.

        Scope is determined by quadrature_rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        Args:
            quadrature_rule: Quadrature rule
            domain: The domain of the integral
            v: the ufl expression
            vaccess: the LNodes expression to access the value in the code
        """
        self.scopes[(domain, quadrature_rule)][v] = vaccess

    def get_var(self, quadrature_rule, domain, v):
        """Lookup ufl expression v in variable scope dicts.

        Scope is determined by quadrature rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        If v is not found in quadrature loop scope, the piecewise
        scope (None) is checked.

        Returns the LNodes expression to access the value in the code.
        """
        if v._ufl_is_literal_:
            return L.ufl_to_lnodes(v)

        # quadrature loop scope
        f = self.scopes[(domain, quadrature_rule)].get(v)

        # piecewise scope
        if f is None:
            f = self.scopes[(None, None)].get(v)
        return f

    def new_temp_symbol(self, basename, dtype=L.DataType.SCALAR):
        """Create a new code symbol named basename + running counter."""
        name = f"{basename}{self.symbol_counters[basename]:d}"
        self.symbol_counters[basename] += 1
        return L.Symbol(name, dtype=dtype)

    def get_temp_symbol(self, tempname, key):
        """Get a temporary symbol."""
        key = (tempname,) + key
        s = self.temp_symbols.get(key)
        defined = s is not None
        if not defined:
            s = self.new_temp_symbol(tempname)
            self.temp_symbols[key] = s
        return s, defined

    def generate(self, domain: basix.CellType):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a
        context that matches a suitable version of the UFC
        tabulate_tensor signatures.
        """
        # Assert that scopes are empty: expecting this to be called only
        # once
        assert not any(d for d in self.scopes.values())

        parts: list[L.LNode] = []

        # Generate the tables of quadrature points and weights
        parts += self.generate_quadrature_tables(domain, self.ir.expression)

        # Generate the tables of basis function values and
        # pre-integrated blocks
        parts += self.generate_element_tables(domain)

        # Generate the tables of geometry data that are needed
        parts += geometry.static_tables(
            self.ir.expression.entity_type,
            self.ir.expression.integrand.values(),
            self.ir.expression.integration_domain_coordinate_element,
        )

        # Loop generation code will produce parts to go before
        # quadloops, to define the quadloops, and to go after the
        # quadloops
        all_preparts = []
        all_quadparts = []

        # Pre-definitions are collected across all quadrature loops to
        # improve re-use and avoid name clashes
        for cell, rule in self.ir.expression.integrand.keys():
            if domain == cell:
                # Generate code to compute piecewise constant scalar factors
                all_preparts += self.generate_piecewise_partition(rule, cell)

                # Generate code to integrate reusable blocks of final
                # element tensor
                all_quadparts += self.generate_quadrature_loop(rule, cell)

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts

        # Drop every unread scalar declaration now the kernel body is whole.
        # Done on the LNode tree, before formatting, so it respects scopes
        # and is shared by all backends.
        parts = prune_dead_scalars(parts, scalar_declaration_names(parts))

        return L.StatementList(parts)

    def generate_quadrature_tables(
        self, domain: basix.CellType, expression: CommonExpressionIR
    ) -> list[L.LNode]:
        """Generate static tables of quadrature points and weights."""
        parts: list[L.LNode] = []
        # No quadrature tables for custom (given argument)
        skip = ufl.custom_integral_types
        if expression.integral_type in skip:
            return parts

        # Loop over quadrature rules. A sum-factorised contraction folds the
        # weights into its 1D tables.
        for (cell, quadrature_rule), integrand in expression.integrand.items():
            if domain == cell and integrand["sum_factorisation"] is None:
                # Generate quadrature weights array
                wsym = self.backend.symbols.weights_table(quadrature_rule)
                parts += [L.ArrayDecl(wsym, values=quadrature_rule.weights, const=True)]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, "Quadrature rules")
        return parts

    def generate_element_tables(self, domain: basix.CellType) -> list[L.LNode]:
        """Generate static tables.

        With precomputed element basis function values in quadrature points.
        """
        parts: list[L.LNode] = []
        tables = self.ir.expression.unique_tables[domain]
        table_types = self.ir.expression.unique_table_types[domain]
        if self.ir.expression.integral_type in ufl.custom_integral_types:
            # Define only piecewise tables
            table_names = [name for name in sorted(tables) if table_types[name] in piecewise_ttypes]
        else:
            # Define all tables
            table_names = sorted(tables)

        for name in table_names:
            table = tables[name]
            parts += self.declare_table(name, table)

        # Add leading comment if there are any tables
        parts = L.commented_code_list(
            parts,
            [
                "Precomputed values of basis functions and precomputations",
                "FE* dimensions: [permutation][entities][points][dofs]",
            ],
        )
        return parts

    def declare_table(self, name, table) -> list[L.ArrayDecl]:
        """Declare a table.

        If the dof dimensions of the table have dof rotations, apply
        these rotations.

        """
        table_symbol = L.Symbol(name, dtype=L.DataType.REAL)
        self.backend.symbols.element_tables[name] = table_symbol
        return [L.ArrayDecl(table_symbol, values=table, const=True)]

    def generate_quadrature_loop(self, quadrature_rule: QuadratureRule, domain: basix.CellType):
        """Generate quadrature loop with for this quadrature_rule."""
        # Generate varying partition
        definitions, intermediates_0 = self.generate_varying_partition(quadrature_rule, domain)

        plan = self.ir.expression.integrand[(domain, quadrature_rule)]["sum_factorisation"]
        if plan is not None:
            return self._sum_factorised_quadrature_loop(
                definitions, intermediates_0, quadrature_rule, domain, plan
            )

        # Generate dofblock parts, some of this will be placed before or after quadloop
        tensor_comp, intermediates_fw = self.generate_dofblock_partition(quadrature_rule, domain)
        assert all([isinstance(tc, L.Section) for tc in tensor_comp])

        # Check if we only have Section objects
        inputs = []
        for definition in definitions:
            assert isinstance(definition, L.Section)
            inputs += definition.output

        # Create intermediates section
        output = []
        declarations = []
        for fw in intermediates_fw:
            assert isinstance(fw, L.VariableDecl)
            output += [fw.symbol]
            # No initialiser: it's unconditionally assigned below, so one
            # would be a dead store. Declared here rather than at that
            # assignment because fw must stay visible to the "Tensor
            # Computation" section that follows, outside this block's scope.
            declarations += [L.VariableDecl(fw.symbol)]
            intermediates_0 += [L.Assign(fw.symbol, fw.value)]
        intermediates = [L.Section("Intermediates", intermediates_0, declarations, inputs, output)]

        iq_symbol = self.backend.symbols.quadrature_loop_index
        iq = create_quadrature_index(quadrature_rule, iq_symbol)

        code = definitions + intermediates + tensor_comp
        code = optimize(code, quadrature_rule)

        split = self._split_quadrature_loop(code, iq, intermediates_fw, quadrature_rule)
        if split is not None:
            return split

        return [L.create_nested_for_loops([iq], code)]

    def _blocks(
        self, quadrature_rule: QuadratureRule, domain: basix.CellType
    ) -> list[tuple[tuple, BlockDataT]]:
        """Return the element tensor blocks of a quadrature rule, in a fixed order."""
        block_contributions = self.ir.expression.integrand[(domain, quadrature_rule)][
            "block_contributions"
        ]
        return [
            (blockmap, blockdata)
            for blockmap, contributions in sorted(block_contributions.items())
            for blockdata in contributions
        ]

    def _sum_factorised_quadrature_loop(
        self,
        definitions: list[L.LNode],
        intermediates: list[L.LNode],
        quadrature_rule: QuadratureRule,
        domain: basix.CellType,
        plan: SumFactorisationPlan,
    ) -> list[L.LNode]:
        """Evaluate the integrand on the tensor product rule, then contract one direction at a time.

        The fused loop updates every entry of the element tensor at every
        quadrature point, costing O(Q^d N^(rd)) for rank r and dimension d.
        With the integrand factors f cached on the whole tensor product
        rule, each direction can be contracted in turn: in 2D, for a
        bilinear form with 1D weights w0 and w1,

            g[iq1] = sum_iq0 f[iq0, iq1] w0[iq0] u0[iq0, i0] v0[iq0, j0]
            A[i0 i1, j0 j1] += sum_iq1 g[iq1] w1[iq1] u1[iq1, i1] v1[iq1, j1]

        which costs O(N^2 Q^2 + N^4 Q). Products of argument tables and the
        1D weights are precomputed in the IR. Definitions annotated with
        `factorize` are evaluated on the whole rule before the integrand.
        """
        integrand = self.ir.expression.integrand[(domain, quadrature_rule)]
        F = integrand["factorization"]
        iq = create_quadrature_index(quadrature_rule, self.backend.symbols.quadrature_loop_index)
        rule_id = quadrature_rule.id()

        # Cache every integrand factor on the whole rule
        caches: dict[int, L.Symbol] = {}
        store: list[L.LNode] = []
        for group in plan.groups:
            for term in group.terms:
                if term.factor not in caches:
                    f = self.get_var(quadrature_rule, domain, F.nodes[term.factor]["expression"])
                    cache = L.Symbol(f"fq{len(caches)}_{rule_id}", dtype=L.DataType.SCALAR)
                    caches[term.factor] = cache
                    store.append(L.Assign(cache[iq.global_index], f))

        prior: list[L.LNode] = []
        pointwise: list[L.Section] = []
        for definition in definitions:
            assert isinstance(definition, L.Section)
            if L.Annotation.factorize in definition.annotations:
                prior += definition.statements
            else:
                pointwise.append(definition)
        inputs = [symbol for definition in pointwise for symbol in definition.output]
        evaluation: list[L.LNode] = [
            *pointwise,
            L.Section("Intermediates", intermediates + store, [], inputs, list(caches.values())),
        ]
        evaluation = optimize(evaluation, quadrature_rule)

        contraction: list[L.LNode] = []
        for group_id, group in enumerate(plan.groups):
            contraction += self._sum_factorised_contraction(plan, group, iq, caches, group_id)

        caches_declarations: list[L.LNode] = [
            L.ArrayDecl(cache, sizes=math.prod(iq.sizes)) for cache in caches.values()
        ]
        return (
            caches_declarations
            + prior
            + [L.create_nested_for_loops([iq], evaluation)]
            + contraction
        )

    def _sum_factorised_contraction(
        self,
        plan: SumFactorisationPlan,
        group: ContractionGroup,
        iq: L.MultiIndex,
        caches: dict[int, L.Symbol],
        group_id: int,
    ) -> list[L.LNode]:
        """Contract cached integrand factors with the argument 1D tables, direction by direction.

        All terms share the same dof loops. Stage arrays shared by several
        terms are computed once. In the last direction, each element tensor
        entry is a sum over the points of the terms' contributions; terms
        sharing the 1D tables of all but the last argument, such as the
        cross terms of a Laplacian on a non-affine cell, are first combined
        over the last argument's table,

            h[q, j] = sum_t g_t[q] B_t[q, j]
            A[i, j] += sum_q Bw[q, i] h[q, j]

        so each entry takes one product per group of terms instead of two
        per term. Blocks receiving the same contributions, such as the
        components of a vector Laplacian, share one sum.
        """
        B_indices = [
            create_dof_index(td, self.backend.symbols.argument_loop_index(i))
            for i, td in enumerate(group.terms[0].tables)
        ]
        dim = iq.dim
        outer_indices, inner_index = B_indices[:-1], B_indices[-1]
        last = dim - 1
        num_points = iq.sizes[last]
        num_inner = inner_index.sizes[last]
        q, j = iq.local_index(last), inner_index.local_index(last)

        def dof_loops(direction: int, body: list[L.LNode], indices) -> list[L.LNode]:
            # Last argument innermost, so the innermost loop is contiguous in A
            for B_index in reversed(indices):
                body = [
                    L.ForRange(B_index.local_index(direction), 0, B_index.sizes[direction], body)
                ]
            return body

        def weighted(direction, tables, quadrature_index, dof_indices):
            """Access the product of 1D tables and the 1D quadrature weights."""
            table = self.backend.symbols.element_tables[plan.weighted_table(direction, tables)]
            return table[(quadrature_index, *dof_indices)]

        def factor(td, direction, quadrature_index, dof_index):
            """Access a 1D table."""
            access, _ = self.backend.access.tensor_factor_access(
                td, direction, quadrature_index, dof_index
            )
            return access

        def total(values: list[L.LExpr]) -> L.LExpr:
            return values[0] if len(values) == 1 else L.Sum(values)

        # level_code[s] runs inside the dof loops of directions <= s. A
        # stage array g_s holds a cache contracted over directions < s, on
        # the points of directions >= s.
        level_code: list[list[L.LNode]] = [[] for _ in range(dim)]
        declarations: list[L.LNode] = []
        stages: dict[tuple, L.Symbol] = {}
        for term in group.terms:
            stage = caches[term.factor]
            keys = term.stage_keys()
            for s in range(last):
                if keys[s + 1] not in stages:
                    stride = math.prod(iq.sizes[s + 1 :])
                    next_stage = L.Symbol(
                        f"{caches[term.factor].name}_{group_id}_g{s + 1}_{len(stages)}",
                        dtype=L.DataType.SCALAR,
                    )
                    stages[keys[s + 1]] = next_stage
                    declarations.append(L.ArrayDecl(next_stage, sizes=stride))
                    m = L.Symbol(f"iqm{s}", dtype=L.DataType.INT)
                    qs = iq.local_index(s)
                    scale = self.new_temp_symbol("basis", dtype=L.DataType.REAL)
                    basis = weighted(
                        s, term.tables, qs, [index.local_index(s) for index in B_indices]
                    )
                    # Accumulate independent entries along the contiguous last
                    # axis, preserving the quadrature order of each entry.
                    width = iq.sizes[-1]
                    indices = [m * width + lane for lane in range(width)]
                    level_code[s] += [
                        L.ForRange(m, 0, stride, [L.Assign(next_stage[m], 0.0)]),
                        L.ForRange(
                            qs,
                            0,
                            iq.sizes[s],
                            [
                                L.VariableDecl(scale, basis),
                                L.ForRange(
                                    m,
                                    0,
                                    stride // width,
                                    [
                                        L.AssignAdd(
                                            next_stage[index], stage[qs * stride + index] * scale
                                        )
                                        for index in indices
                                    ],
                                ),
                            ],
                        ),
                    ]
                stage = stages[keys[s + 1]]

        def last_stage(term) -> L.Symbol:
            """Return the cache of a term contracted over all but the last direction."""
            return stages[term.stage_keys()[last]] if last > 0 else caches[term.factor]

        # pre_rows runs inside the dof loops of directions < last, before
        # the last direction's; row_code inside those of the other arguments
        pre_rows: list[L.LNode] = []
        row_code: list[L.LNode] = []
        A = self.backend.symbols.element_tensor
        for row in group.rows:
            single = [(last_stage(term), term.tables) for term in row.single]
            combined = []
            for members in row.combined:
                tables = members[0].tables
                outer = [
                    weighted(last, tables[:1], q, [outer_indices[0].local_index(last)]),
                    *[
                        factor(td, last, q, index.local_index(last))
                        for td, index in zip(tables[1:-1], outer_indices[1:])
                    ],
                ]
                # h[q, j] = sum over the group of stage[q] * table[q, j]
                h = self.new_temp_symbol("h")
                declarations.append(L.ArrayDecl(h, sizes=num_points * num_inner))
                value = total(
                    [
                        last_stage(member)[q] * factor(member.tables[-1], last, q, j)
                        for member in members
                    ]
                )
                pre_rows.append(
                    L.ForRange(
                        q,
                        0,
                        num_points,
                        [L.ForRange(j, 0, num_inner, [L.Assign(h[q * num_inner + j], value)])],
                    )
                )
                combined.append((outer, h))

            def update(dof):
                """Return one entry's contribution at the current quadrature point."""
                values = [
                    L.float_product([*outer, h[q * num_inner + dof]]) for outer, h in combined
                ]
                indices = [index.local_index(last) for index in outer_indices] + [dof]
                values += [
                    stage[q] * weighted(last, tables, q, indices) for stage, tables in single
                ]
                return total(values)

            acc = self.new_temp_symbol("acc")
            declarations.append(L.ArrayDecl(acc, sizes=num_inner))
            # Each dof has an independent accumulator. Exposing the fixed-size
            # row lets the compiler pack the updates without reassociating sums.
            row_body: list[L.LNode] = [
                L.ForRange(j, 0, num_inner, [L.Assign(acc[j], L.LiteralFloat(0.0))]),
                L.ForRange(
                    q,
                    0,
                    num_points,
                    [L.AssignAdd(acc[dof], update(L.LiteralInt(dof))) for dof in range(num_inner)],
                ),
            ]
            store = []
            for target in row.targets:
                A_indices = self._element_tensor_indices(
                    target.blockmap, target.blockdata, B_indices
                )
                A_index = L.MultiIndex(A_indices, self.ir.expression.tensor_shape)
                store.append(L.AssignAdd(A[A_index], acc[j]))
            row_body.append(L.ForRange(j, 0, num_inner, store))
            row_code += row_body

        body = pre_rows + dof_loops(last, row_code, outer_indices)
        for s in reversed(range(last)):
            body = dof_loops(s, level_code[s] + body, B_indices)
        return declarations + body

    def _element_tensor_indices(
        self, blockmap: tuple, blockdata: BlockDataT, B_indices: list[L.MultiIndex]
    ) -> list[L.LExpr]:
        """Return the indices into the element tensor of a block's argument dofs."""
        A_indices = []
        for i, index in enumerate(B_indices):
            tabledata = blockdata.ma_data[i].tabledata
            offset = tabledata.offset
            if len(blockmap[i]) == 1:
                A_indices.append(index.global_index + offset)
            else:
                A_indices.append(tabledata.block_size * index.global_index + offset)
        return A_indices

    def _split_quadrature_loop(
        self,
        code: list[L.LNode],
        iq: L.MultiIndex,
        intermediates_fw: list[L.VariableDecl],
        quadrature_rule: QuadratureRule,
    ):
        """Split the quadrature loop into evaluation and contraction loops.

        The fused loop interleaves evaluation of the varying quantities
        (which carries any math functions of the spatial coordinate) with
        the tensor contraction, whose strided table reads stop the compiler
        vectorising the loop -- leaving calls to `sin` and friends scalar.
        Caching the `fw` values and contracting in a second loop leaves the
        evaluation loop free of those reads, so it vectorises onto the SIMD
        math library.

        The `fw` cache is tiled, not sized to the full point count: a form
        with many `fw` values and a large rule could otherwise put megabytes
        on the stack. Tiling bounds stack use to a fixed budget and keeps
        each tile resident in L1d across both passes.

        Returns `None` if the split does not apply, in which case the
        caller emits a single fused loop.
        """
        if quadrature_rule is None or not intermediates_fw:
            return None

        # Only single-index (non tensor-product) rules for now
        if quadrature_rule.has_tensor_factors or iq.dim != 1:
            return None

        evaluation = [c for c in code if getattr(c, "name", None) != "Tensor Computation"]
        contraction = [c for c in code if getattr(c, "name", None) == "Tensor Computation"]
        if not contraction or not evaluation:
            return None

        # The split only pays for itself when there's a transcendental math
        # call to vectorise -- most integrands (mass matrices, plain
        # Poisson, etc.) have none, and for those it's a net loss.
        if not _contains_vectorizable_mathfunction(evaluation):
            return None
        # On x86-64 Linux the evaluation loop can call libmvec. Elsewhere,
        # only split when the contraction is expensive enough to amortise
        # the fw cache without that benefit.
        if not _supports_vector_math_loop_split() and (
            _contraction_work(contraction) < _UNSUPPORTED_VECTOR_MATH_MIN_CONTRACTION_WORK
        ):
            return None

        num_points = quadrature_rule.weights.size
        index = iq.local_index(0)

        # Budget the whole fw cache at 32 KiB, comfortably inside a typical
        # L1d, using a conservative 16 bytes/scalar (covers complex128).
        n_fw = len(intermediates_fw)
        tile = _fw_cache_tile_size(n_fw, num_points)
        if tile is None:
            return None

        caches = {
            fw.symbol.name: L.Symbol(f"{fw.symbol.name}_q", fw.symbol.dtype)
            for fw in intermediates_fw
        }
        declarations: list[L.LNode] = [L.ArrayDecl(cache, sizes=tile) for cache in caches.values()]

        if tile >= num_points:
            # Whole rule fits in one tile: no need for an outer tiling loop.
            store = [
                L.Assign(L.ArrayAccess(caches[fw.symbol.name], [index]), fw.symbol)
                for fw in intermediates_fw
            ]
            load = [
                L.VariableDecl(fw.symbol, L.ArrayAccess(caches[fw.symbol.name], [index]))
                for fw in intermediates_fw
            ]
            return declarations + [
                L.create_nested_for_loops([iq], evaluation + store),
                L.create_nested_for_loops([iq], load + contraction),
            ]

        num_tiles = -(-num_points // tile)  # ceil division
        tile_index = L.Symbol(f"{index.name}_tile", L.DataType.INT)
        tile_base = L.Symbol(f"{index.name}_base", L.DataType.INT)
        tile_end = L.Symbol(f"{index.name}_end", L.DataType.INT)
        local_offset = L.Sub(index, tile_base)

        tile_end_expr = L.Conditional(
            L.LT(L.Add(tile_base, L.LiteralInt(tile)), L.LiteralInt(num_points)),
            L.Add(tile_base, L.LiteralInt(tile)),
            L.LiteralInt(num_points),
        )
        setup: list[L.LNode] = [
            L.VariableDecl(tile_base, L.Mul(tile_index, L.LiteralInt(tile))),
            L.VariableDecl(tile_end, tile_end_expr),
        ]

        store = [
            L.Assign(L.ArrayAccess(caches[fw.symbol.name], [local_offset]), fw.symbol)
            for fw in intermediates_fw
        ]
        load = [
            L.VariableDecl(fw.symbol, L.ArrayAccess(caches[fw.symbol.name], [local_offset]))
            for fw in intermediates_fw
        ]

        tile_body = setup + [
            L.ForRange(index, tile_base, tile_end, evaluation + store),
            L.ForRange(index, tile_base, tile_end, load + contraction),
        ]
        return declarations + [L.ForRange(tile_index, 0, num_tiles, tile_body)]

    def generate_piecewise_partition(self, quadrature_rule, domain: basix.CellType):
        """Generate a piecewise partition."""
        # Get annotated graph of factorisation
        F = self.ir.expression.integrand[(domain, quadrature_rule)]["factorization"]
        arraysymbol = L.Symbol(f"sp_{quadrature_rule.id()}", dtype=L.DataType.SCALAR)
        return self.generate_partition(arraysymbol, F, "piecewise", None, None)

    def generate_varying_partition(self, quadrature_rule, domain: basix.CellType):
        """Generate a varying partition."""
        # Get annotated graph of factorisation
        F = self.ir.expression.integrand[(domain, quadrature_rule)]["factorization"]
        arraysymbol = L.Symbol(f"sv_{quadrature_rule.id()}", dtype=L.DataType.SCALAR)
        return self.generate_partition(arraysymbol, F, "varying", quadrature_rule, domain)

    def generate_partition(self, symbol, F, mode, quadrature_rule, domain):
        """Generate a partition."""
        definitions = []
        intermediates = []

        for i, attr in F.nodes.items():
            if attr["status"] != mode:
                continue
            v = attr["expression"]

            # Generate code only if the expression is not already in cache
            if not self.get_var(quadrature_rule, domain, v):
                if v._ufl_is_literal_:
                    vaccess = L.ufl_to_lnodes(v)
                elif mt := attr.get("mt"):
                    tabledata = attr.get("tr")

                    # Backend specific modified terminal translation
                    vaccess = self.backend.access.get(mt, tabledata, quadrature_rule)
                    vdef = self.backend.definitions.get(mt, tabledata, quadrature_rule, vaccess)

                    # A sum-factorised definition is a list of sections
                    for section in vdef if isinstance(vdef, list) else [vdef]:
                        assert isinstance(section, L.Section)
                        # Only add if definition is unique.
                        # This can happen when using sub-meshes
                        if section not in definitions:
                            definitions += [section]
                else:
                    # Get previously visited operands
                    vops = [self.get_var(quadrature_rule, domain, op) for op in v.ufl_operands]
                    dtype = extract_dtype(v, vops)

                    # Mapping UFL operator to target language
                    self._ufl_names.add(v._ufl_handler_name_)
                    vexpr = L.ufl_to_lnodes(v, *vops)

                    j = len(intermediates)
                    vaccess = L.Symbol(f"{symbol.name}_{j}", dtype=dtype)
                    intermediates.append(L.VariableDecl(vaccess, vexpr))

                # Store access node for future reference
                self.set_var(quadrature_rule, domain, v, vaccess)

        # Optimize definitions
        definitions = optimize(definitions, quadrature_rule)

        # Fold scalar chains that round-trip through exact powers of two
        # (e.g. a symmetric gradient's 1/2 undoing an earlier *2) into a
        # plain alias. Must run before reciprocal-CSE, which would otherwise
        # hide the power-of-two divisors this looks for behind a symbol.
        intermediates = power_of_two_cse_statements(intermediates)

        # Collapse repeated divisions by the same divisor (e.g. an affine
        # cell's pseudo-inverse Jacobian dividing every entry by the same
        # determinant) into one reciprocal and a multiply each.
        intermediates = reciprocal_cse_statements(intermediates, name_prefix=f"recip_{symbol.name}")

        # Record each piecewise two-factor product so generate_block_parts
        # can fuse a shared factor across block entries with the quadrature
        # weight once (see _piecewise_mul_operands). Read after the CSE
        # passes so a reciprocal-CSE-rewritten division counts too.
        if mode == "piecewise":
            for stmt in intermediates:
                if isinstance(stmt, L.VariableDecl) and isinstance(stmt.value, L.Mul):
                    self._piecewise_mul_operands[stmt.symbol.name] = (
                        stmt.value.lhs,
                        stmt.value.rhs,
                    )

        return definitions, intermediates

    def generate_dofblock_partition(
        self,
        quadrature_rule: QuadratureRule,
        domain: basix.CellType,
    ):
        """Generate a dofblock partition."""
        quadparts = []
        blocks = self._blocks(quadrature_rule, domain)

        block_groups = collections.defaultdict(list)

        # Group loops by blockmap, in Vector elements each component has
        # a different blockmap
        for blockmap, blockdata in blocks:
            scalar_blockmap = []
            assert len(blockdata.ma_data) == len(blockmap)
            for i, b in enumerate(blockmap):
                bs = blockdata.ma_data[i].tabledata.block_size
                offset = blockdata.ma_data[i].tabledata.offset
                # Only sum-factorisation tensor factors lack these, and they
                # are never the table of a block's modified argument
                assert bs is not None and offset is not None
                b = tuple([(idx - offset) // bs for idx in b])
                scalar_blockmap.append(b)
            block_groups[tuple(scalar_blockmap)].append(blockdata)

        intermediates = []
        for blockmap in block_groups:
            block_quadparts, intermediate = self.generate_block_parts(
                quadrature_rule,
                domain,
                blockmap,
                block_groups[blockmap],
            )
            intermediates += intermediate

            # Add computations
            quadparts.extend(block_quadparts)

        return quadparts, intermediates

    def get_arg_factors(self, blockdata, block_rank, quadrature_rule, domain, iq, indices):
        """Get arg factors."""
        arg_factors = []
        tables = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            scope = self.ir.expression.integrand[(domain, quadrature_rule)]["modified_arguments"]
            mt = scope[mad.ma_index]
            arg_tables = []

            # Translate modified terminal to code
            # TODO: Move element table access out of backend?
            #       Not using self.backend.access.argument() here
            #       now because it assumes too much about indices.

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = 1
            else:
                # Assuming B sparsity follows element table sparsity
                arg_factor, arg_tables = self.backend.access.table_access(
                    td, self.ir.expression.entity_type, mt.restriction, iq, indices[i]
                )

            tables += arg_tables
            arg_factors.append(arg_factor)

        return arg_factors, tables

    def generate_block_parts(
        self,
        quadrature_rule: QuadratureRule,
        domain: basix.CellType,
        blockmap: tuple,
        blocklist: list[BlockDataT],
    ):
        """Generate and return code parts for a given block.

        Returns parts occurring before, inside, and after the quadrature
        loop identified by the quadrature rule.

        Should be called with quadrature_rule=None for
        quadloop-independent blocks.
        """
        # The parts to return
        quadparts: list[L.LNode] = []
        intermediates: list[L.LNode] = []
        tables = []
        vars = []

        # RHS expressions grouped by LHS "dofmap"
        rhs_expressions = collections.defaultdict(list)

        block_rank = len(blockmap)
        iq_symbol = self.backend.symbols.quadrature_loop_index
        iq = create_quadrature_index(quadrature_rule, iq_symbol)

        A_shape = self.ir.expression.tensor_shape

        # Detect block entries whose piecewise values share one common
        # factor (e.g. several metric-tensor entries scaled by the same
        # |det J|). Multiplying that factor into the quadrature weight once,
        # rather than into every already-scaled entry, is cheaper whenever
        # there are more sharing entries than quadrature points. This is a
        # genuine re-association GCC can't perform itself without
        # -ffast-math, since it changes rounding.
        #
        # Skip it when every argument table in this block is itself
        # quadrature-point-independent (piecewise_ttypes) and there is more
        # than one quadrature point: the per-entry values are then fully
        # loop-invariant, and GCC/clang hoist the unfused multiply-by-weight
        # out of the loop entirely -- an indirection this fusion adds (a
        # shared factor computed fresh each iteration) can stop the compiler
        # seeing that the whole loop is redundant, costing far more than the
        # fusion saves (see PR #865's hyperelasticity_residual_p1 finding).
        # With a single quadrature point there's no loop to hoist away, so
        # the fusion's reduced op count is a clean win regardless.
        fused_factor: dict[int, tuple[L.LExpr, L.LExpr]] = {}
        if (
            quadrature_rule is not None
            and self.ir.expression.integral_type not in ufl.custom_integral_types
            and len(blocklist) > quadrature_rule.weights.size
            and (
                quadrature_rule.weights.size == 1
                or any(
                    ttype not in piecewise_ttypes
                    for blockdata in blocklist
                    for ttype in blockdata.ttypes
                )
            )
        ):
            F = self.ir.expression.integrand[(domain, quadrature_rule)]["factorization"]
            factor_indices = []
            decomps: list[tuple[L.LExpr, L.LExpr]] = []
            for blockdata in blocklist:
                if len(blockdata.factor_indices_comp_indices) > 1:
                    decomps = []
                    break
                factor_index = blockdata.factor_indices_comp_indices[0][0]
                v = F.nodes[factor_index]["expression"]
                f = self.get_var(quadrature_rule, domain, v)
                decomp = (
                    self._piecewise_mul_operands.get(f.name) if isinstance(f, L.Symbol) else None
                )
                if decomp is None:
                    decomps = []
                    break
                factor_indices.append(factor_index)
                decomps.append(decomp)

            if decomps:
                rhs_keys = {expr_key(d[1]) for d in decomps}
                lhs_keys = {expr_key(d[0]) for d in decomps}
                shared, bases = None, None
                if len(rhs_keys) == 1:
                    shared, bases = decomps[0][1], [d[0] for d in decomps]
                elif len(lhs_keys) == 1:
                    shared, bases = decomps[0][0], [d[1] for d in decomps]
                if shared is not None and bases is not None:
                    fused_factor = dict(zip(factor_indices, zip(bases, [shared] * len(bases))))

        for blockdata in blocklist:
            B_indices = []
            for i in range(block_rank):
                table_ref = blockdata.ma_data[i].tabledata
                symbol = self.backend.symbols.argument_loop_index(i)
                index = create_dof_index(table_ref, symbol)
                B_indices.append(index)

            if self.ir.part == TensorPart.diagonal and block_rank == 2:
                assert len(A_shape) == 1
                B_indices = [B_indices[0], B_indices[0]]

            ttypes = blockdata.ttypes
            if "zeros" in ttypes:
                raise RuntimeError(
                    "Not expecting zero arguments to be left in dofblock generation."
                )

            if len(blockdata.factor_indices_comp_indices) > 1:
                raise RuntimeError("Code generation for non-scalar integrals unsupported")

            # We have scalar integrand here, take just the factor index
            factor_index = blockdata.factor_indices_comp_indices[0][0]

            # Get factor expression
            F = self.ir.expression.integrand[(domain, quadrature_rule)]["factorization"]

            v = F.nodes[factor_index]["expression"]
            f = self.get_var(quadrature_rule, domain, v)

            # Quadrature weight was removed in representation, add it back now
            if self.ir.expression.integral_type in ufl.custom_integral_types:
                weights = self.backend.symbols.custom_weights_table
                weight = weights[iq.global_index]
            else:
                weights = self.backend.symbols.weights_table(quadrature_rule)
                weight = weights[iq.global_index]

            if factor_index in fused_factor:
                base, shared = fused_factor[factor_index]
                scale_key = (quadrature_rule, expr_key(shared))
                combined, defined = self.get_temp_symbol("fwscale", scale_key)
                if not defined:
                    intermediates += [L.VariableDecl(combined, L.float_product([shared, weight]))]
                f, weight = base, combined

            # Define fw = f * weight
            fw_rhs = L.float_product([f, weight])
            if not isinstance(fw_rhs, L.Product):
                fw = fw_rhs
            else:
                # Define and cache scalar temp variable
                key = (quadrature_rule, factor_index, blockdata.all_factors_piecewise)
                fw, defined = self.get_temp_symbol("fw", key)
                if not defined:
                    input = [f, weight]
                    # filter only L.Symbol in input
                    input = [i for i in input if isinstance(i, L.Symbol)]
                    output = [fw]

                    # assert input and output are Symbol objects
                    assert all(isinstance(i, L.Symbol) for i in input)
                    assert all(isinstance(o, L.Symbol) for o in output)

                    intermediates += [L.VariableDecl(fw, fw_rhs)]

            var = fw if isinstance(fw, L.Symbol) else fw.array
            vars += [var]
            assert not blockdata.transposed, "Not handled yet"

            # Fetch code to access modified arguments
            arg_factors, table = self.get_arg_factors(
                blockdata, block_rank, quadrature_rule, domain, iq, B_indices
            )
            tables += table
            # Define B_rhs = fw * arg_factors
            insert_rank = block_rank
            if self.ir.part == TensorPart.diagonal:
                insert_rank = 1
                B_indices = [B_indices[0]]
            B_rhs = L.float_product([fw] + arg_factors)

            A_indices = self._element_tensor_indices(blockmap, blockdata, B_indices[:insert_rank])
            rhs_expressions[tuple(A_indices)].append(B_rhs)

        # List of statements to keep in the inner loop
        keep = collections.defaultdict(list)

        for indices in rhs_expressions:
            keep[indices] = rhs_expressions[indices]

        body: list[L.LNode] = []

        A = self.backend.symbols.element_tensor
        for indices in keep:
            multi_index = L.MultiIndex(list(indices), A_shape)
            for expression in keep[indices]:
                body.append(L.AssignAdd(A[multi_index], expression))

        # Nest with the last tensor index innermost when an entry sums more
        # than one term (e.g. gradient-gradient) -- GCC's vectoriser wants
        # this. A single-term entry (e.g. a plain mass matrix) prefers the
        # opposite, row-innermost order instead.
        multi_term = any(len(v) > 1 for v in keep.values())
        nest_indices = B_indices if multi_term else B_indices[::-1]
        body = [L.create_nested_for_loops(nest_indices, body)]
        input = [*vars, *tables]
        output = [A]

        # Make sure we don't have repeated symbols in input
        input = list(set(input))

        # assert input and output are Symbol objects
        assert all(isinstance(i, L.Symbol) for i in input)
        assert all(isinstance(o, L.Symbol) for o in output)

        annotations = []
        if len(B_indices) > 1:
            annotations.append(L.Annotation.licm)

        quadparts += [L.Section("Tensor Computation", body, [], input, output, annotations)]

        return quadparts, intermediates
