# Copyright (C) 2015-2024 Martin Sandve Alnæs, Michal Habera, Igor Baratta, Chris Richardson
#
# Modified by Jørgen S. Dokken, 2024, 2026
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Integral generator."""

import collections
import itertools
import logging
import math
import platform
import sys
from numbers import Integral
from typing import Any, cast

import basix
import numpy as np
import ufl

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration import geometry
from ffcx.codegeneration.definitions import create_dof_index, create_quadrature_index
from ffcx.codegeneration.optimizer import (
    expr_key,
    factor_scalar_sums,
    optimize,
    power_of_two_cse_statements,
    prune_dead_scalars,
    prune_unused_tables,
    reciprocal_cse_statements,
    referenced_symbol_names,
    scalar_declaration_names,
)
from ffcx.ir.elementtables import piecewise_ttypes
from ffcx.ir.integral import BlockDataT, CommonExpressionIR, TensorPart
from ffcx.ir.representation import IntegralIR
from ffcx.ir.representationutils import QuadratureRule

logger = logging.getLogger("ffcx")


_FW_CACHE_BUDGET_BYTES = 32768
_MAX_LOCAL_TENSOR_BYTES = 65536
_LOCAL_TENSOR_ALIGNMENT_BYTES = 32
_MAX_PREINTEGRATED_TABLE_BYTES = 65536
_MIN_PARTIAL_PREINTEGRATION_WORK = 64
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
    if num_points <= max_tile:
        return num_points

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


def _children(node) -> list:
    """Return the direct sub-nodes of a statement or expression."""
    if isinstance(node, list):
        return node
    if isinstance(node, L.Section):
        return [*node.declarations, *node.statements]
    if isinstance(node, L.StatementList):
        return list(node.statements)
    if isinstance(node, L.ForRange):
        return list(node.body.statements)
    if isinstance(node, L.VariableDecl):
        return [] if node.value is None else [node.value]
    if isinstance(node, (L.ArrayDecl, L.Comment)):
        return []
    if isinstance(node, L.Statement):
        return [node.expr]
    if isinstance(node, L.ArrayAccess):
        return [node.array, *node.indices]
    if isinstance(node, L.MultiIndex):
        return list(node.symbols)
    if isinstance(node, L.BinOp):
        return [node.lhs, node.rhs]
    if isinstance(node, (L.NaryOp, L.MathFunction)):
        return list(node.args)
    if isinstance(node, L.PrefixUnaryOp):
        return [node.arg]
    if isinstance(node, L.Conditional):
        return [node.condition, node.true, node.false]
    return []


def _flat_statements(statements) -> list[L.LNode]:
    """Flatten nested statement lists into one list of statements."""
    out: list[L.LNode] = []
    for s in statements:
        if isinstance(s, L.StatementList):
            out += _flat_statements(s.statements)
        else:
            out.append(s)
    return out


def _references(node, symbol: L.Symbol) -> bool:
    """Return whether *node* refers to *symbol* anywhere."""
    if isinstance(node, L.Symbol):
        return node.name == symbol.name
    return any(_references(child, symbol) for child in _children(node))


def _writes_symbol(node, symbol: L.Symbol) -> bool:
    """Return whether *node* assigns to an entry of the array *symbol*."""
    if isinstance(node, L.AssignOp):
        return isinstance(node.lhs, L.ArrayAccess) and node.lhs.array == symbol
    return any(_writes_symbol(child, symbol) for child in _children(node))


def _is_index(node, index: L.Symbol) -> bool:
    """Return whether *node* is the loop index itself.

    Accepts the index wrapped as a one-term sum, or as a product with unit
    literals, as flattened multi-indices produce.
    """
    if isinstance(node, L.MultiIndex):
        return node.dim == 1 and _is_index(node.symbols[0], index)
    if isinstance(node, (L.Sum, L.Product)):
        unit = isinstance(node, L.Product)
        args = [a for a in node.args if not (unit and isinstance(a, L.LiteralInt) and a.value == 1)]
        return len(args) == 1 and _is_index(args[0], index)
    if isinstance(node, L.Mul):
        for a, b in ((node.lhs, node.rhs), (node.rhs, node.lhs)):
            if isinstance(a, L.LiteralInt) and a.value == 1:
                return _is_index(b, index)
        return False
    return isinstance(node, L.Symbol) and node.name == index.name


def _tables_indexed_last_by(
    exprs: list[L.LExpr], index: L.Symbol, tables: dict[str, L.Symbol]
) -> set[str] | None:
    """Return element tables read with *index* as their final subscript.

    ``None`` means *index* is also used another way, so padding the loop
    over *index* past the table extent would not be safe.
    """
    found: set[str] = set()

    def visit(node) -> bool:
        if isinstance(node, L.Symbol):
            return node.name != index.name
        if isinstance(node, L.ArrayAccess) and _references(node, index):
            name = node.array.name
            if (
                name not in tables
                or not _is_index(node.indices[-1], index)
                or any(_references(i, index) for i in node.indices[:-1])
            ):
                return False
            found.add(name)
            return True
        return all(visit(child) for child in _children(node))

    return found if all(visit(e) for e in exprs) else None


def _map_expressions(node, fn):
    """Rebuild a statement or expression tree, replacing nodes with ``fn``.

    ``fn`` returns a replacement for a node or ``None`` to keep it and
    rebuild its operands.
    """
    if isinstance(node, list):
        return [_map_expressions(n, fn) for n in node]
    replacement = fn(node)
    if replacement is not None:
        return replacement
    if isinstance(node, L.Section):
        return L.Section(
            node.name,
            _map_expressions(list(node.statements), fn),
            _map_expressions(list(node.declarations), fn),
            node.input,
            node.output,
            node.annotations,
        )
    if isinstance(node, L.StatementList):
        return L.StatementList(_map_expressions(list(node.statements), fn))
    if isinstance(node, L.ForRange):
        return L.ForRange(
            node.index, node.begin, node.end, _map_expressions(list(node.body.statements), fn)
        )
    if isinstance(node, L.VariableDecl):
        value = None if node.value is None else _map_expressions(node.value, fn)
        return L.VariableDecl(node.symbol, value)
    if isinstance(node, (L.ArrayDecl, L.Comment)):
        return node
    if isinstance(node, L.Statement):
        return L.Statement(_map_expressions(node.expr, fn))
    if isinstance(node, L.ArrayAccess):
        return L.ArrayAccess(node.array, [_map_expressions(i, fn) for i in node.indices])
    if isinstance(node, L.MathFunction):
        return L.MathFunction(node.function, [_map_expressions(a, fn) for a in node.args])
    if isinstance(node, L.NaryOp):
        return type(node)([_map_expressions(a, fn) for a in node.args])
    if isinstance(node, L.BinOp):
        return type(node)(_map_expressions(node.lhs, fn), _map_expressions(node.rhs, fn))
    if isinstance(node, L.PrefixUnaryOp):
        return type(node)(_map_expressions(node.arg, fn))
    if isinstance(node, L.Conditional):
        return L.Conditional(
            _map_expressions(node.condition, fn),
            _map_expressions(node.true, fn),
            _map_expressions(node.false, fn),
        )
    return node


def _contraction_nest(section: L.Section, A: L.Symbol, arguments: list[L.Symbol]):
    """Return the loop nest of a contraction section that updates A.

    Returns ``(position, indices, ranges, updates)`` for a section whose only
    update of A is one perfect nest of one or two zero-based loops with
    ``A[...] += expr`` statements, with ``indices`` ordered (test, trial).
    Returns ``None`` for any other form.
    """
    nests = [
        (n, s)
        for n, s in enumerate(section.statements)
        if isinstance(s, L.ForRange) and _writes_symbol(s, A)
    ]
    if len(nests) != 1 or any(
        _writes_symbol(s, A) for s in section.statements if not isinstance(s, L.ForRange)
    ):
        return None
    position, nest = nests[0]
    loops = []
    node: L.LNode = nest
    while isinstance(node, L.ForRange):
        loops.append(node)
        statements = _flat_statements(node.body.statements)
        if len(statements) == 1 and isinstance(statements[0], L.ForRange):
            node = statements[0]
        else:
            break
    inner = _flat_statements(loops[-1].body.statements)
    if (
        len(loops) not in (1, 2)
        or not all(isinstance(s, L.Statement) and isinstance(s.expr, L.AssignAdd) for s in inner)
        or not all(isinstance(loop.index, L.Symbol) for loop in loops)
        or not all(isinstance(loop.begin, L.LiteralInt) and loop.begin.value == 0 for loop in loops)
        or not all(isinstance(loop.end, L.LiteralInt) for loop in loops)
    ):
        return None
    ranges = {loop.index: int(loop.end.value) for loop in loops}
    if len(loops) == 1:
        indices = [loops[0].index]
    elif set(ranges) == set(arguments):
        indices = list(arguments)
    else:
        return None
    updates = [s.expr for s in inner if isinstance(s, L.Statement)]
    if not all(
        isinstance(u.lhs, L.ArrayAccess) and u.lhs.array == A and not _references(u.rhs, A)
        for u in updates
    ):
        return None
    return position, indices, ranges, updates


def _evaluate_index(node, values: dict[str, int]) -> int:
    """Evaluate an integer index expression for given loop index values."""
    if isinstance(node, (list, tuple)):
        if len(node) != 1:
            raise ValueError("Expecting a flattened index.")
        return _evaluate_index(node[0], values)
    if isinstance(node, L.LiteralInt):
        return int(node.value)
    if isinstance(node, L.Symbol):
        return values[node.name]
    if isinstance(node, L.MultiIndex):
        return _evaluate_index(node.global_index, values)
    if isinstance(node, L.Sum):
        return sum(_evaluate_index(a, values) for a in node.args)
    if isinstance(node, L.Product):
        return math.prod(_evaluate_index(a, values) for a in node.args)
    if isinstance(node, L.Add):
        return _evaluate_index(node.lhs, values) + _evaluate_index(node.rhs, values)
    if isinstance(node, L.Sub):
        return _evaluate_index(node.lhs, values) - _evaluate_index(node.rhs, values)
    if isinstance(node, L.Mul):
        return _evaluate_index(node.lhs, values) * _evaluate_index(node.rhs, values)
    if isinstance(node, L.Neg):
        return -_evaluate_index(node.arg, values)
    raise ValueError(f"Cannot evaluate index {node!r}.")


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


def _arithmetic_work(node) -> int:
    """Estimate scalar arithmetic in one evaluation of a generated section."""
    if isinstance(node, list):
        return sum(_arithmetic_work(child) for child in node)
    if isinstance(node, L.Section):
        return _arithmetic_work(node.declarations) + _arithmetic_work(node.statements)
    if isinstance(node, L.StatementList):
        return _arithmetic_work(node.statements)
    if isinstance(node, L.ForRange):
        if not isinstance(node.begin, L.LiteralInt) or not isinstance(node.end, L.LiteralInt):
            return 0
        return int(max(0, node.end.value - node.begin.value)) * _arithmetic_work(node.body)
    if isinstance(node, L.VariableDecl):
        return _arithmetic_work(node.value)
    if isinstance(node, (L.ArrayDecl, L.Comment)):
        return 0
    if isinstance(node, L.Statement):
        return _arithmetic_work(node.expr)
    if isinstance(node, L.Assign):
        return _arithmetic_work(node.rhs)
    if isinstance(node, L.BinOp):
        return 1 + _arithmetic_work(node.lhs) + _arithmetic_work(node.rhs)
    if isinstance(node, L.NaryOp):
        return max(0, len(node.args) - 1) + sum(_arithmetic_work(arg) for arg in node.args)
    if isinstance(node, L.MathFunction):
        return 1 + sum(_arithmetic_work(arg) for arg in node.args)
    if isinstance(node, L.PrefixUnaryOp):
        return _arithmetic_work(node.arg)
    if isinstance(node, L.Conditional):
        return _arithmetic_work(node.condition) + max(
            _arithmetic_work(node.true), _arithmetic_work(node.false)
        )
    return 0


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

        # Scalar symbol name -> (lhs, rhs) of its defining two-factor
        # product, if any. Lets generate_block_parts spot several block
        # entries sharing one factor (e.g. several metric-tensor entries all
        # scaled by |det J|) and fuse it with the quadrature weight once.
        self._scalar_mul_operands: dict[str, tuple[L.LExpr, L.LExpr]] = {}

        # Set of counters used for assigning names to intermediate
        # variables
        self.symbol_counters: dict[str, int] = collections.defaultdict(int)

        # Element table name -> padded length of its final (dof) axis, for
        # contiguous loops that read zero entries past the last dof.
        self._padded_tables: dict[str, int] = {}

        # Element table name -> padded point count of a copy with point and
        # dof axes exchanged, read by loops that vectorise across points.
        # Weights table name -> padded point count.
        self._transposed_tables: dict[str, int] = {}
        self._padded_weights: dict[str, int] = {}

        # Scalar symbol name -> id of its value, equal for symbols with
        # structurally identical definitions (see _value_key). Used to prove
        # that two factors are identical.
        self._value_ids: dict[str, int] = {}
        self._interned: dict[tuple, int] = {}
        self._polynomials: dict[str, dict[tuple[int, ...], Any]] = {}

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

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
        name = f"{basename}{self.symbol_counters[basename]:d}"
        self.symbol_counters[basename] += 1
        return L.Symbol(name, dtype=L.DataType.SCALAR)

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

        parts: list[L.LNode] = [
            L.Comment("Accumulate this integral into A; existing tensor entries are incremented.")
        ]

        # Generate the tables of quadrature points and weights
        parts += self.generate_quadrature_tables(domain, self.ir.expression)

        # Generate the tables of basis function values and
        # pre-integrated blocks
        element_tables = self.generate_element_tables(domain)
        parts += element_tables

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

        # Zero-pad the dof axis of tables read by padded contiguous loops.
        # Existing reads index real dofs only, so their values are unchanged.
        for decl in element_tables:
            if isinstance(decl, L.ArrayDecl) and decl.symbol.name in self._padded_tables:
                size = self._padded_tables[decl.symbol.name]
                values = np.asarray(decl.values)
                pad = [(0, 0)] * (values.ndim - 1) + [(0, size - values.shape[-1])]
                decl.values = np.pad(values, pad)
                decl.sizes = decl.values.shape

        # Point-major copies of tables for loops across quadrature points.
        # Padded points repeat the last point, so every value stays finite.
        transposed: list[L.LNode] = []
        unique_tables = cast(dict[str, np.ndarray], self.ir.expression.unique_tables[domain])
        for name, points in sorted(self._transposed_tables.items()):
            values = np.swapaxes(np.asarray(unique_tables[name]), -1, -2)
            pad = [(0, 0)] * (values.ndim - 1) + [(0, points - values.shape[-1])]
            symbol = L.Symbol(f"{name}_T", dtype=L.DataType.REAL)
            transposed.append(
                L.ArrayDecl(symbol, values=np.pad(values, pad, mode="edge"), const=True)
            )
        if transposed:
            position = len(parts)
            if element_tables:
                position = next(k for k, p in enumerate(parts) if p is element_tables[-1]) + 1
            parts[position:position] = L.commented_code_list(
                transposed, "Point-major copies of tables for evaluation across points."
            )
        for decl in parts:
            if isinstance(decl, L.ArrayDecl) and decl.symbol.name in self._padded_weights:
                values = np.asarray(decl.values)
                size = self._padded_weights[decl.symbol.name]
                decl.values = np.pad(values, (0, size - values.shape[-1]))
                decl.sizes = decl.values.shape

        # Drop every unread scalar declaration now the kernel body is whole.
        # Done on the LNode tree, before formatting, so it respects scopes
        # and is shared by all backends.
        parts = prune_dead_scalars(parts, scalar_declaration_names(parts))
        parts = prune_unused_tables(parts)

        return L.StatementList(parts)

    def generate_quadrature_tables(self, domain: basix.CellType, expression: CommonExpressionIR):
        """Generate static tables of quadrature points and weights."""
        parts: list[L.LNode] = []
        # No quadrature tables for custom (given argument)
        skip = ufl.custom_integral_types
        if expression.integral_type in skip:
            return parts

        # Loop over quadrature rules
        for (cell, quadrature_rule), _ in expression.integrand.items():
            if domain == cell and not self._quadrature_invariant(quadrature_rule, domain):
                # Generate quadrature weights array
                wsym = self.backend.symbols.weights_table(quadrature_rule)
                parts += [L.ArrayDecl(wsym, values=quadrature_rule.weights, const=True)]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(
            parts,
            [
                "Quadrature weights on the reference integration entity.",
                "Geometry scaling is applied when assembling the tensor.",
            ],
        )
        return parts

    def generate_element_tables(self, domain: basix.CellType):
        """Generate static tables.

        With precomputed element basis function values in quadrature points.
        """
        parts = []
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
                "Precomputed reference basis values and derivatives.",
                "FE* dimensions: [permutation][entities][points][dofs]",
            ],
        )
        return parts

    def declare_table(self, name, table):
        """Declare a table.

        If the dof dimensions of the table have dof rotations, apply
        these rotations.

        """
        table_symbol = L.Symbol(name, dtype=L.DataType.REAL)
        self.backend.symbols.element_tables[name] = table_symbol
        return [L.ArrayDecl(table_symbol, values=table, const=True)]

    def _generate_preintegrated_tensor(
        self, quadrature_rule: QuadratureRule, domain: basix.CellType
    ) -> tuple[list[L.LNode], set[int]] | None:
        """Precompute reference contractions for blocks with point-independent factors.

        Each table holds an argument table or product of two tables summed over
        quadrature points. Runtime assembly contracts these reference tensors
        with the cellwise factors, avoiding repeated output updates for each
        point. Contributions sharing a factor and output map are added into one
        reference table. The table budget bounds data for large argument spaces.
        Return the generated code and identities of the integrated IR blocks so
        the remaining quadrature contraction can omit those exact contributions.
        """
        if (
            self.ir.expression.integral_type != "cell"
            or self.ir.part != TensorPart.full
            or quadrature_rule.has_tensor_factors
            or quadrature_rule.weights.size < 2
            or self._quadrature_invariant(quadrature_rule, domain)
        ):
            return None
        integrand = self.ir.expression.integrand[(domain, quadrature_rule)]
        blocks = [
            b
            for contributions in integrand["block_contributions"].values()
            for b in contributions
            if len(b.ma_data) in (1, 2)
            and b.all_factors_piecewise
            and not b.transposed
            and len(b.factor_indices_comp_indices) == 1
            and all(
                not ma.tabledata.has_tensor_factorisation
                and not ma.tabledata.is_permuted
                and ma.tabledata.values.shape[:2] == (1, 1)
                for ma in b.ma_data
            )
        ]
        if not blocks:
            return None
        # A partial contraction adds another pass over the output. Keep small
        # mixed forms with varying argument tables fused unless enough terms
        # are removed to pay for that pass and the reference-table reads.
        # Point-independent argument tables can always collapse their loop.
        num_blocks = sum(
            len(contributions) for contributions in integrand["block_contributions"].values()
        )
        if (
            len(blocks) < num_blocks
            and len(blocks) * (quadrature_rule.weights.size - 1) < _MIN_PARTIAL_PREINTEGRATION_WORK
            and any(not ma.tabledata.is_piecewise for block in blocks for ma in block.ma_data)
        ):
            return None
        pairs = {
            tuple(ma.tabledata.name for ma in b.ma_data): tuple(ma.tabledata for ma in b.ma_data)
            for b in blocks
        }
        table_bytes = sum(
            math.prod(t.values.shape[-1] for t in pair) * 8 for pair in pairs.values()
        )
        if table_bytes > _MAX_PREINTEGRATED_TABLE_BYTES:
            return None

        tables: dict[tuple[str, ...], np.ndarray] = {}
        for names, pair in pairs.items():
            argument_values = [
                np.broadcast_to(
                    table.values[0, 0], (len(quadrature_rule.weights), table.values.shape[-1])
                )
                for table in pair
            ]
            contraction = "q,qi->i" if len(pair) == 1 else "q,qi,qj->ij"
            tables[names] = np.einsum(contraction, quadrature_rule.weights, *argument_values)

        groups: dict[tuple, dict[tuple, tuple[L.LExpr, np.ndarray]]] = collections.defaultdict(dict)
        term_counts: collections.Counter[tuple] = collections.Counter()
        F = integrand["factorization"]
        A = self.backend.symbols.element_tensor
        shape = self.ir.expression.tensor_shape
        indices = [self.backend.symbols.argument_loop_index(k) for k in range(len(shape))]
        for block in blocks:
            data = [ma.tabledata for ma in block.ma_data]
            assert all(t.block_size is not None and t.offset is not None for t in data)
            factor = self.get_var(
                quadrature_rule,
                domain,
                F.nodes[block.factor_indices_comp_indices[0][0]]["expression"],
            )
            sizes = tuple(t.values.shape[-1] for t in data)
            output = tuple(t.block_size * index + t.offset for t, index in zip(data, indices))
            values = tables[tuple(t.name for t in data)]
            # A shared cellwise factor can multiply the sum of reference tables.
            # Merge only contributions with the same output map and dimensions.
            entries = groups[sizes, output]
            term_counts[sizes, output] += 1
            key = expr_key(factor)
            # Preserve the unrolling of small constant-argument updates when
            # they are followed by a remaining quadrature contraction.
            if len(blocks) < num_blocks and all(t.is_piecewise for t in data):
                key = (key, term_counts[sizes, output])
            if key in entries:
                previous_factor, previous_values = entries[key]
                entries[key] = (previous_factor, previous_values + values)
            else:
                entries[key] = (factor, values)

        declarations: list[L.LNode] = []
        symbols: dict[tuple, L.Symbol] = {}
        loops: dict[tuple[int, ...], list[L.LNode]] = collections.defaultdict(list)
        table_bytes = 0
        for (sizes, output), entries in groups.items():
            # A full contiguous tensor has no row boundaries to preserve. One
            # linear loop lets vectorizers use contiguous loads and stores.
            # Keep existing loops for unmerged contractions: flattening can
            # worsen unrolling for small mass tensors or many table streams.
            flat = (
                sizes == tuple(shape)
                and output == tuple(indices)
                and 1 < len(entries) < term_counts[sizes, output]
            )
            table_indices = [indices[0]] if flat else indices
            terms = []
            for factor, values in entries.values():
                if flat:
                    values = values.ravel()
                key = (values.shape, values.tobytes())
                if key not in symbols:
                    symbol = L.Symbol(f"RI_{quadrature_rule.id()}_{len(symbols)}", L.DataType.REAL)
                    symbols[key] = symbol
                    declarations.append(L.ArrayDecl(symbol, values=values, const=True))
                    table_bytes += values.nbytes
                terms.append(factor * L.ArrayAccess(symbols[key], table_indices))
            target = A[indices[0]] if flat else A[L.MultiIndex(list(output), shape)]
            loop_sizes = (math.prod(sizes),) if flat else sizes
            loops[loop_sizes].append(L.AssignAdd(target, L.Sum(terms)))
        if table_bytes > _MAX_PREINTEGRATED_TABLE_BYTES:
            return None
        code: list[L.LNode] = []
        for sizes, body in loops.items():
            for index, size in reversed(list(zip(indices, sizes))):
                body = [L.ForRange(index, 0, size, body)]
            if len(sizes) == 1 and len(shape) == 2:
                comment = "Traverse the full tensor in contiguous storage order."
            else:
                comment = "Map argument dof indices to element tensor entries."
            code.extend([L.Comment(comment), *body])
        declarations = L.commented_code_list(
            declarations,
            [
                "Reference integrals (RI_*): argument basis values or products integrated",
                "with quadrature weights at code generation time.",
                "Contributions sharing a scalar factor and output map may share a summed table.",
            ],
        )
        code = L.commented_code_list(
            code,
            "Contract reference integrals with factors constant over quadrature points; add to A.",
        )
        return declarations + code, {id(block) for block in blocks}

    def generate_quadrature_loop(self, quadrature_rule: QuadratureRule, domain: basix.CellType):
        """Generate quadrature loop with for this quadrature_rule."""
        preintegrated = self._generate_preintegrated_tensor(quadrature_rule, domain)
        preintegrated_code, preintegrated_blocks = preintegrated or ([], set())
        blocks = self.ir.expression.integrand[(domain, quadrature_rule)]["block_contributions"]
        if preintegrated_blocks and all(
            id(block) in preintegrated_blocks
            for contributions in blocks.values()
            for block in contributions
        ):
            return preintegrated_code

        # Generate varying partition
        definitions, intermediates_0 = self.generate_varying_partition(quadrature_rule, domain)

        # Generate dofblock parts, some of this will be placed before or after quadloop
        tensor_comp, intermediates_fw = self.generate_dofblock_partition(
            quadrature_rule, domain, preintegrated_blocks
        )
        assert all([isinstance(tc, L.Section) for tc in tensor_comp])

        invariant = self._quadrature_invariant(quadrature_rule, domain)
        prologue: list[L.LNode] = []
        epilogue: list[L.LNode] = []
        if not invariant:
            tensor_comp, prologue, epilogue = self._localize_tensor_updates(tensor_comp)
            # Blocks recovered by symmetry no longer read their factors.
            needed = referenced_symbol_names(list(tensor_comp))
            kept = []
            for fw in reversed(intermediates_fw):
                if fw.symbol.name in needed:
                    kept.append(fw)
                    needed |= referenced_symbol_names([fw])
            intermediates_fw = kept[::-1]

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

        if invariant:
            return preintegrated_code + L.commented_code_list(
                code,
                [
                    "The integrand is constant across quadrature points.",
                    "Apply the summed quadrature weight once, without a point loop.",
                ],
            )

        split = self._split_quadrature_loop(code, iq, intermediates_fw, quadrature_rule)
        if split is not None:
            return preintegrated_code + prologue + split + epilogue

        return (
            preintegrated_code
            + prologue
            + [
                L.Comment(
                    f"Quadrature rule {quadrature_rule.id()}: evaluate and accumulate at "
                    f"{quadrature_rule.weights.size} points."
                ),
                L.create_nested_for_loops([iq], code),
            ]
            + epilogue
        )

    def _localize_tensor_updates(
        self, sections: list[L.Section]
    ) -> tuple[list[L.Section], list[L.LNode], list[L.LNode]]:
        """Accumulate quadrature-loop tensor updates in local arrays.

        For contractions updating several blocks of A, such as the
        components of blocked elements, updates of each block accumulate into
        their own zero-initialised array, with the trial (or only) argument
        index innermost and contiguous. Blocked elements otherwise update
        interleaved entries of A at a stride, which defeats vectorisation, and
        every quadrature point would read and write A again. After the
        quadrature loop, one pass adds the local arrays into A. Performing
        the arithmetic before touching A also lets it overlap with any memory
        traffic for A. The inner loop runs to a 32-byte multiple when every
        read indexed by it is a zero-padded element table.

        A block is skipped if its terms are exactly those of another block,
        or those with test and trial tables exchanged while its entries are
        the transposed entries of that block, with identical scalar factors
        (see _equivalent_blocks). Its update of A then reads the other local
        array, transposed if required.

        Returns the transformed sections, declarations to place before the
        quadrature loop and the final update of A to place after it.
        Sections with fewer than three distinct blocks or without a plain
        contraction loop nest are kept, and all sections are kept if the
        local arrays exceed a fixed stack budget.
        """
        A = self.backend.symbols.element_tensor
        arguments = [self.backend.symbols.argument_loop_index(k) for k in range(2)]
        scalar_type = np.dtype(self.backend.access.options["scalar_type"])
        scalar_bytes = scalar_type.itemsize
        if np.issubdtype(scalar_type, np.complexfloating):
            real_bytes = scalar_bytes // 2
        else:
            real_bytes = scalar_bytes
        width = _LOCAL_TENSOR_ALIGNMENT_BYTES // real_bytes

        out: list[L.Section] = []
        prologue: list[L.LNode] = []
        epilogue: list[L.LNode] = []
        padded: dict[str, int] = {}
        total_bytes = 0
        for section in sections:
            nest = _contraction_nest(section, A, arguments)
            if nest is None:
                out.append(section)
                continue
            position, indices, ranges, updates = nest

            # Group the updates by element tensor entry map
            groups: dict[tuple, tuple[L.ArrayAccess, list[L.LExpr]]] = {}
            for update in updates:
                groups.setdefault(expr_key(update.lhs), (update.lhs, []))[1].append(update.rhs)
            equivalent = self._equivalent_blocks(list(groups.values()), indices, ranges)
            if len(groups) - len(equivalent) < 3:
                # A single block already updates contiguous rows of A, in a
                # loop order chosen for it. A few blocks leave too little work
                # per entry to pay for the extra pass over local arrays, and
                # small loop bodies can lead compilers to vectorise across
                # rows with gathers and scatters.
                out.append(section)
                continue

            innermost = indices[-1]
            tables = _tables_indexed_last_by(
                [u.rhs for u in updates], innermost, self.backend.symbols.element_tables
            )
            size = ranges[innermost]
            if tables is not None and size % width:
                size += width - size % width
                for name in tables:
                    padded[name] = max(padded.get(name, 0), size)
            sizes = [ranges[index] for index in indices[:-1]] + [size]
            total_bytes += (len(groups) - len(equivalent)) * math.prod(sizes) * scalar_bytes
            if total_bytes > _MAX_LOCAL_TENSOR_BYTES:
                return sections, [], []

            body: list[L.LNode] = []
            scatter: list[L.LNode] = []
            locals_: dict[int, L.Symbol] = {}
            for g, (lhs, rhs) in enumerate(groups.values()):
                if g in equivalent:
                    continue
                local = self.new_temp_symbol("A_local_")
                locals_[g] = local
                # A single initial value, nested to the array rank, zero-fills
                # the whole array.
                prologue.append(L.ArrayDecl(local, sizes=sizes, values=np.zeros((1,) * len(sizes))))
                body += [L.AssignAdd(L.ArrayAccess(local, indices), r) for r in rhs]
            for g, (lhs, _) in enumerate(groups.values()):
                source, transpose = equivalent.get(g, (g, False))
                access = L.ArrayAccess(locals_[source], indices[::-1] if transpose else indices)
                scatter.append(L.AssignAdd(lhs, access))
            for index, end in reversed(list(zip(indices, sizes))):
                body = [L.ForRange(index, 0, end, body)]
            for index in reversed(indices):
                scatter = [L.ForRange(index, 0, ranges[index], scatter)]
            epilogue += scatter
            statements = list(section.statements)
            statements[position] = body[0]
            out.append(
                L.Section(
                    section.name,
                    statements,
                    section.declarations,
                    section.input,
                    section.output,
                    section.annotations,
                )
            )

        if not prologue:
            return sections, [], []
        for name, size in padded.items():
            self._padded_tables[name] = max(self._padded_tables.get(name, 0), size)
        prologue = L.commented_code_list(
            prologue,
            "Local element tensor blocks, accumulated over quadrature points.",
        )
        epilogue = L.commented_code_list(
            epilogue, "Add the local blocks into A in a single pass after the quadrature loop."
        )
        return out, prologue, epilogue

    def _record_value(self, symbol: L.Symbol, value: L.LExpr, expand: bool = True) -> None:
        """Record the definition of a scalar symbol for value comparisons.

        With ``expand``, small polynomial definitions are expanded over their
        atoms, which include unexpanded symbols, so equal polynomials
        evaluated in different orders share an id. Coefficients must match
        exactly. Definitions are recorded in order, so operands already have
        ids and keys never recurse through chains of definitions.
        """
        polynomial = self._polynomial(value) if expand else None
        if polynomial is None:
            self._value_ids[symbol.name] = self._value_key(value)
        elif list(polynomial.values()) == [1] and len(monomial := next(iter(polynomial))) == 1:
            # A single atom, such as an alias, keeps the id of that atom.
            self._value_ids[symbol.name] = monomial[0]
        else:
            self._polynomials[symbol.name] = polynomial
            terms = sorted((monomial, repr(c)) for monomial, c in polynomial.items())
            self._value_ids[symbol.name] = self._intern(("polynomial", *terms))

    def _polynomial(self, node) -> dict[tuple[int, ...], Any] | None:
        """Expand an expression into a polynomial over atom ids, if small.

        Returns ``{monomial: coefficient}``, with sorted atom ids as the
        monomial, or ``None`` when the expansion exceeds a fixed size.
        """
        max_terms = 64

        def add(a, b, sign=1):
            out = dict(a)
            for k, v in b.items():
                out[k] = out.get(k, 0) + sign * v
            return {k: v for k, v in out.items() if v != 0}

        def mul(a, b):
            if len(a) * len(b) > max_terms:
                return None
            out: dict[tuple[int, ...], Any] = {}
            for ka, va in a.items():
                for kb, vb in b.items():
                    k = tuple(sorted(ka + kb))
                    out[k] = out.get(k, 0) + va * vb
            return {k: v for k, v in out.items() if v != 0}

        def expand(node):
            if isinstance(node, (L.LiteralFloat, L.LiteralInt)):
                return {(): node.value} if node.value != 0 else {}
            if isinstance(node, L.Symbol) and node.name in self._polynomials:
                return self._polynomials[node.name]
            if isinstance(node, (L.Add, L.Sub, L.Sum)):
                args = node.args if isinstance(node, L.Sum) else [node.lhs, node.rhs]
                out: dict | None = {}
                for k, arg in enumerate(args):
                    term = expand(arg)
                    if term is None or out is None:
                        return None
                    out = add(out, term, -1 if isinstance(node, L.Sub) and k else 1)
                return out if len(out) <= max_terms else None
            if isinstance(node, L.Neg):
                term = expand(node.arg)
                return None if term is None else {k: -v for k, v in term.items()}
            if isinstance(node, (L.Mul, L.Product)):
                args = node.args if isinstance(node, L.Product) else [node.lhs, node.rhs]
                out = {(): 1}
                for arg in args:
                    term = expand(arg)
                    if term is None or out is None:
                        return None
                    out = mul(out, term)
                return out
            if isinstance(node, L.Div) and isinstance(node.rhs, (L.LiteralFloat, L.LiteralInt)):
                term = expand(node.lhs)
                return None if term is None else {k: v / node.rhs.value for k, v in term.items()}
            return {(self._value_key(node),): 1}

        result: dict[tuple[int, ...], Any] | None = expand(node)
        return result

    def _value_key(self, node) -> int:
        """Return an id shared by structurally identical expressions.

        Recorded scalar symbols take the id of their definition, and the
        operands of commutative operations are ordered. Ids are interned so
        large expression graphs are keyed in linear time.
        """
        if isinstance(node, L.Symbol):
            if node.name not in self._value_ids:
                self._value_ids[node.name] = self._intern(("symbol", node.name))
            return self._value_ids[node.name]
        if isinstance(node, (L.LiteralFloat, L.LiteralInt)):
            return self._intern(("literal", node.value))
        children = [self._value_key(c) for c in _children(node)]
        if isinstance(node, (L.Sum, L.Product, L.Add, L.Mul)):
            children.sort()
        return self._intern((type(node).__name__, getattr(node, "function", None), *children))

    def _intern(self, key: tuple) -> int:
        """Return a unique integer for a hashable key."""
        return self._interned.setdefault(key, len(self._interned))

    def _equivalent_blocks(
        self,
        groups: list[tuple[L.ArrayAccess, list[L.LExpr]]],
        indices: list[L.Symbol],
        ranges: dict[L.Symbol, int],
    ) -> dict[int, tuple[int, bool]]:
        """Map blocks equal to an earlier block, or its transpose, to that block.

        Values are compared through their terms, each a product of scalar
        factors, one test table and one trial table, so two blocks are equal
        when their terms are the same and transposed when their tables are
        exchanged. A transposed block also needs a square tensor, with entry
        maps evaluated for every dof pair to check that they are transposed.
        Returns ``{block: (earlier block, transposed)}``.
        """
        if len(indices) != 2 or self.ir.part != TensorPart.full:
            return {}
        test, trial = indices
        tables = self.backend.symbols.element_tables

        def table_key(node, index):
            if (
                isinstance(node, L.ArrayAccess)
                and node.array.name in tables
                and _references(node.indices[-1], index)
                and not any(_references(i, index) for i in node.indices[:-1])
            ):
                return (node.array.name, *(self._value_key(i) for i in node.indices[:-1]))
            return None

        def signature(rhs: list[L.LExpr]):
            terms: collections.Counter[tuple] = collections.Counter()
            for term in rhs:
                if not isinstance(term, L.Product):
                    return None
                rows, cols, scalars = [], [], []
                for arg in term.args:
                    if _references(arg, test):
                        rows.append(table_key(arg, test))
                    elif _references(arg, trial):
                        cols.append(table_key(arg, trial))
                    else:
                        scalars.append(self._value_key(arg))
                if len(rows) != 1 or len(cols) != 1 or None in rows + cols:
                    return None
                terms[(tuple(sorted(scalars)), rows[0], cols[0])] += 1
            return terms

        signatures = [signature(rhs) for _, rhs in groups]
        if any(s is None for s in signatures):
            return {}
        equivalent: dict[int, tuple[int, bool]] = {}
        for h, sh in enumerate(signatures):
            for g in range(h):
                if g not in equivalent and signatures[g] == sh:
                    equivalent[h] = (g, False)
                    break

        shape = self.ir.expression.tensor_shape
        n = ranges[test]
        if len(shape) != 2 or shape[0] != shape[1] or ranges[trial] != n or n > 64:
            return equivalent

        def entries(lhs):
            try:
                return tuple(
                    divmod(_evaluate_index(lhs.indices, {test.name: i, trial.name: j}), shape[1])
                    for i in range(n)
                    for j in range(n)
                )
            except (KeyError, TypeError, ValueError):
                return None

        maps = [entries(lhs) for lhs, _ in groups]
        if any(m is None for m in maps):
            return equivalent
        # Entries of a block's transpose, in the same (test, trial) dof order
        transposed: dict[tuple, int] = {}
        for g, m in enumerate(maps):
            assert m is not None
            t = tuple((m[j * n + i][1], m[j * n + i][0]) for i in range(n) for j in range(n))
            transposed.setdefault(t, g)
        for h, m in enumerate(maps):
            source = transposed.get(m)
            if source is None or source >= h or source in equivalent or h in equivalent:
                continue
            sg, sh = signatures[source], signatures[h]
            assert sg is not None and sh is not None
            swapped = collections.Counter({(s, c, r): k for (s, r, c), k in sg.items()})
            if swapped == sh:
                equivalent[h] = (source, True)
        return equivalent

    def _quadrature_invariant(
        self, quadrature_rule: QuadratureRule, domain: basix.CellType
    ) -> bool:
        """Check whether every block is independent of the quadrature point."""
        if self.ir.expression.integral_type in ufl.custom_integral_types:
            return False
        integrand = self.ir.expression.integrand[(domain, quadrature_rule)]
        blocks = [b for values in integrand["block_contributions"].values() for b in values]
        return bool(blocks) and all(
            b.all_factors_piecewise and all(t in piecewise_ttypes for t in b.ttypes) for b in blocks
        )

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
        math library. For sufficiently large unrolled geometry evaluations,
        the same separation permits SIMD arithmetic across quadrature points
        without interleaved tensor updates.

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
        contraction_inputs = referenced_symbol_names(contraction)
        intermediates_fw = [fw for fw in intermediates_fw if fw.symbol.name in contraction_inputs]
        if not intermediates_fw:
            return None

        if not _contains_vectorizable_mathfunction(evaluation):
            # Separate unrolled geometry evaluation so SIMD can run across
            # points. Tiny rules, retained reductions and contractions that
            # dominate evaluation do not amortise the additional cache pass.
            jacobians = [
                section for section in evaluation if getattr(section, "name", None) == "Jacobian"
            ]
            reductions = [
                section
                for section in evaluation
                if getattr(section, "name", None) in ("Jacobian", "Coefficient")
            ]
            if (
                not jacobians
                or quadrature_rule.weights.size < 8
                or len(intermediates_fw) > 12
                or any(_contraction_work(section) for section in reductions)
                or _arithmetic_work(evaluation) < _contraction_work(contraction)
                # A cheap coefficient evaluation followed by a tiny load
                # does not reliably amortise the extra cache pass.
                or (
                    any(getattr(section, "name", None) == "Coefficient" for section in reductions)
                    and _arithmetic_work(evaluation) < 64
                    and _contraction_work(contraction) < 16
                )
            ):
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

        evaluation, paddable, transposed = self._point_major(evaluation, index, quadrature_rule)
        points = num_points
        scalar_type = np.dtype(self.backend.access.options["scalar_type"])
        real_bytes = scalar_type.itemsize
        if np.issubdtype(scalar_type, np.complexfloating):
            real_bytes //= 2
        width = _LOCAL_TENSOR_ALIGNMENT_BYTES // real_bytes
        padded = -(-num_points // width) * width
        # Run whole SIMD vectors across points. Factors at the padded points
        # are cached but never contracted. Limit the padding, since a loop
        # that keeps scalar math calls evaluates every padded point.
        affordable = 4 * (padded - num_points) <= num_points
        if tile >= num_points and paddable and transposed and affordable:
            points = padded
            weights = self.backend.symbols.weights_table(quadrature_rule)
            self._padded_weights[weights.name] = points
        for name in transposed:
            self._transposed_tables[name] = max(self._transposed_tables.get(name, 0), points)

        caches = {
            fw.symbol.name: L.Symbol(f"{fw.symbol.name}_q", fw.symbol.dtype)
            for fw in intermediates_fw
        }
        cache_size = points if tile >= num_points else tile
        declarations: list[L.LNode] = [
            L.ArrayDecl(cache, sizes=cache_size) for cache in caches.values()
        ]

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
                L.Comment("Evaluate and cache weighted scalar factors at quadrature points."),
                L.ForRange(index, 0, points, evaluation + store),
                L.Comment(
                    "Reuse cached factors to contract argument tables and accumulate into A."
                ),
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
            L.Comment("Evaluate weighted scalar factors for the current tile of points."),
            L.ForRange(index, tile_base, tile_end, evaluation + store),
            L.Comment("Reuse the tile cache to contract argument tables and accumulate into A."),
            L.ForRange(index, tile_base, tile_end, load + contraction),
        ]
        return declarations + [
            L.Comment("Tile quadrature points to bound temporary storage for weighted factors."),
            L.ForRange(tile_index, 0, num_tiles, tile_body),
        ]

    def _point_major(
        self, evaluation: list[L.LNode], index: L.Symbol, quadrature_rule: QuadratureRule
    ) -> tuple[list[L.LNode], bool, set[str]]:
        """Read element tables point-major in a loop over quadrature points.

        Table reads ``FE[...][iq][k]`` become ``FE_T[...][k][iq]``, so
        consecutive points are contiguous for SIMD across points. Returns
        the rewritten code, whether every read indexed by the point is a
        table or the weights (so the loop may run past the last point), and
        the transposed table names.
        """
        tables = self.backend.symbols.element_tables
        weights = None
        if self.ir.expression.integral_type not in ufl.custom_integral_types:
            weights = self.backend.symbols.weights_table(quadrature_rule).name
        paddable = True
        transposed: set[str] = set()

        def replace(node):
            nonlocal paddable
            if isinstance(node, L.Symbol) and node.name == index.name:
                paddable = False
                return node
            if not isinstance(node, L.ArrayAccess) or not _references(node, index):
                return None
            name = node.array.name
            *outer, point, dof = node.indices if len(node.indices) >= 2 else (None, None)
            if (
                name in tables
                and point is not None
                and _is_index(point, index)
                and not _references(dof, index)
                and not any(_references(i, index) for i in outer)
            ):
                transposed.add(name)
                table = L.Symbol(f"{name}_T", dtype=L.DataType.REAL)
                return L.ArrayAccess(table, [*outer, dof, point])
            if not (name == weights and len(node.indices) == 1):
                paddable = False
            return node

        return _map_expressions(evaluation, replace), paddable, transposed

    def generate_piecewise_partition(self, quadrature_rule, domain: basix.CellType):
        """Generate a piecewise partition."""
        # Get annotated graph of factorisation
        F = self.ir.expression.integrand[(domain, quadrature_rule)]["factorization"]
        arraysymbol = L.Symbol(f"sp_{quadrature_rule.id()}", dtype=L.DataType.SCALAR)
        definitions, intermediates = self.generate_partition(
            arraysymbol, F, "piecewise", None, None
        )
        intermediates = L.commented_code_list(
            intermediates,
            [
                "Scalar factors (sp_*): geometry transforms, coefficient expressions and other",
                "quantities constant across quadrature points of this integration entity.",
            ],
        )
        return definitions, intermediates

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

                    if vdef:
                        assert isinstance(vdef, L.Section)
                    # Only add if definition is unique.
                    # This can happen when using sub-meshes
                    if vdef not in definitions:
                        definitions += [vdef]
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

        # Record each two-factor product so generate_block_parts
        # can fuse a shared factor across block entries with the quadrature
        # weight once (see _scalar_mul_operands). Read after the CSE
        # passes so a reciprocal-CSE-rewritten division counts too.
        for stmt in intermediates:
            if isinstance(stmt, L.VariableDecl) and isinstance(stmt.value, L.Mul):
                self._scalar_mul_operands[stmt.symbol.name] = (
                    stmt.value.lhs,
                    stmt.value.rhs,
                )

        intermediates = factor_scalar_sums(intermediates)

        # Record scalar definitions so the tensor contraction can recognise
        # factors with identical values.
        # Point-varying values are expanded over the cellwise ones, which stay
        # atoms to bound the expansion.
        for stmt in intermediates:
            if isinstance(stmt, L.VariableDecl) and stmt.value is not None:
                self._record_value(stmt.symbol, stmt.value, expand=mode == "varying")

        return definitions, intermediates

    def generate_dofblock_partition(
        self,
        quadrature_rule: QuadratureRule,
        domain: basix.CellType,
        preintegrated_blocks: set[int],
    ):
        """Generate a dofblock partition, excluding already integrated blocks."""
        block_contributions = self.ir.expression.integrand[(domain, quadrature_rule)][
            "block_contributions"
        ]
        quadparts = []
        blocks = [
            (blockmap, blockdata)
            for blockmap, contributions in sorted(block_contributions.items())
            for blockdata in contributions
            if id(blockdata) not in preintegrated_blocks
        ]

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
        rhs_expressions: dict[tuple[L.LExpr, ...], list[L.LExpr]] = collections.defaultdict(list)

        block_rank = len(blockmap)
        iq_symbol = self.backend.symbols.quadrature_loop_index
        iq = create_quadrature_index(quadrature_rule, iq_symbol)

        A_shape = self.ir.expression.tensor_shape

        unroll = all(
            td.is_piecewise
            and td.is_uniform
            and not td.is_permuted
            and not td.has_tensor_factorisation
            and td.values.shape[-1] <= 4
            for b in blocklist
            for td in (ma.tabledata for ma in b.ma_data)
        )
        quadrature_invariant = self._quadrature_invariant(quadrature_rule, domain)

        # Detect block entries whose scalar values share one common
        # factor (e.g. several metric-tensor entries scaled by the same
        # |det J|). Multiplying that factor into the quadrature weight once,
        # rather than into every already-scaled entry, is cheaper whenever
        # there are more sharing entries than quadrature points. This is a
        # genuine re-association GCC can't perform itself without
        # -ffast-math, since it changes rounding.
        # Varying factors save a multiplication per sharing entry at every
        # point, so two entries suffice regardless of the rule size.
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
            and (
                len(blocklist) > quadrature_rule.weights.size
                or (len(blocklist) > 1 and any(not b.all_factors_piecewise for b in blocklist))
            )
            and (
                quadrature_rule.weights.size == 1
                or any(not b.all_factors_piecewise for b in blocklist)
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
                decomp = self._scalar_mul_operands.get(f.name) if isinstance(f, L.Symbol) else None
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
            elif quadrature_invariant:
                weight = L.LiteralFloat(math.fsum(quadrature_rule.weights))
            else:
                weights = self.backend.symbols.weights_table(quadrature_rule)
                weight = weights[iq.global_index]

            if factor_index in fused_factor:
                base, shared = fused_factor[factor_index]
                scale_key = (quadrature_rule, expr_key(shared))
                combined, defined = self.get_temp_symbol("fwscale", scale_key)
                if not defined:
                    value = L.float_product([shared, weight])
                    self._record_value(combined, value)
                    intermediates += [L.VariableDecl(combined, value)]
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

                    self._record_value(fw, fw_rhs)
                    intermediates += [L.VariableDecl(fw, fw_rhs)]

            var = fw if isinstance(fw, L.Symbol) else fw.array
            vars += [var]
            assert not blockdata.transposed, "Not handled yet"

            if unroll:
                data = [ma.tabledata for ma in blockdata.ma_data]
                ranges = [range(td.values.shape[-1]) for td in data]
                for dofs in itertools.product(*ranges):
                    if self.ir.part == TensorPart.diagonal and len(dofs) == 2:
                        if dofs[0] != dofs[1]:
                            continue
                    values = [td.values[0, 0, 0, i] for td, i in zip(data, dofs)]
                    if any(value == 0.0 for value in values):
                        continue
                    output_indices = []
                    for td, i in zip(data, dofs):
                        assert td.block_size is not None and td.offset is not None
                        output_indices.append(L.LiteralInt(td.block_size * i + td.offset))
                    if self.ir.part == TensorPart.diagonal:
                        output_indices = output_indices[:1]
                    rhs_expressions[tuple(output_indices)].append(L.float_product([fw, *values]))
                continue

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

            A_indices = []
            for i in range(insert_rank):
                index = B_indices[i]
                tabledata = blockdata.ma_data[i].tabledata
                offset = tabledata.offset
                if len(blockmap[i]) == 1:
                    A_indices.append(index.global_index + offset)
                else:
                    block_size = blockdata.ma_data[i].tabledata.block_size
                    A_indices.append(block_size * index.global_index + offset)
            rhs_expressions[tuple(A_indices)].append(B_rhs)

        # List of statements to keep in the inner loop
        keep = collections.defaultdict(list)

        for indices in rhs_expressions:
            keep[indices] = rhs_expressions[indices]

        body: list[L.LNode] = []

        A = self.backend.symbols.element_tensor
        for indices in keep:
            multi_index = L.MultiIndex(list(indices), A_shape)
            if unroll or block_rank == 1:
                body.append(L.AssignAdd(A[multi_index], L.Sum(keep[indices])))
            else:
                for expression in keep[indices]:
                    body.append(L.AssignAdd(A[multi_index], expression))

        # Nest with the last tensor index innermost when an entry sums more
        # than one term (e.g. gradient-gradient) -- GCC's vectoriser wants
        # this. A single-term entry (e.g. a plain mass matrix) prefers the
        # opposite, row-innermost order instead.
        multi_term = any(len(v) > 1 for v in keep.values())
        nest_indices = B_indices if multi_term else B_indices[::-1]
        if not unroll:
            body = [L.create_nested_for_loops(nest_indices, body)]
        input = [*vars, *tables]
        output = [A]

        # Make sure we don't have repeated symbols in input
        input = list(set(input))

        # assert input and output are Symbol objects
        assert all(isinstance(i, L.Symbol) for i in input)
        assert all(isinstance(o, L.Symbol) for o in output)

        annotations = []
        if len(B_indices) > 1 and not unroll:
            annotations.append(L.Annotation.licm)

        quadparts += [L.Section("Tensor Computation", body, [], input, output, annotations)]

        return quadparts, intermediates
