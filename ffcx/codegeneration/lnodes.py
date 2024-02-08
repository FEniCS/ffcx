# Copyright (C) 2013-2023 Martin Sandve Alnæs, Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""LNodes.

LNodes is intended as a minimal generic language description.
Formatting is done later, depending on the target language.

Supported:
 Floating point (and complex) and integer variables and multidimensional arrays
 Range loops
 Simple arithmetic, +-*/
 Math operations
 Logic conditions
 Comments
Not supported:
 Pointers
 Function Calls
 Flow control (if, switch, while)
 Booleans
 Strings
"""

import numbers
from collections.abc import Sequence
from enum import Enum
from typing import Optional

import numpy as np
import ufl


class PRECEDENCE:
    """An enum-like class for operator precedence levels."""

    HIGHEST = 0
    LITERAL = 0
    SYMBOL = 0
    SUBSCRIPT = 2

    NOT = 3
    NEG = 3

    MUL = 4
    DIV = 4

    ADD = 5
    SUB = 5

    LT = 7
    LE = 7
    GT = 7
    GE = 7
    EQ = 8
    NE = 8
    AND = 11
    OR = 12
    CONDITIONAL = 13
    ASSIGN = 13
    LOWEST = 15


def is_zero_lexpr(lexpr):
    """Check if an expression is zero."""
    return (isinstance(lexpr, LiteralFloat) and lexpr.value == 0.0) or (
        isinstance(lexpr, LiteralInt) and lexpr.value == 0
    )


def is_one_lexpr(lexpr):
    """Check if an expression is one."""
    return (isinstance(lexpr, LiteralFloat) and lexpr.value == 1.0) or (
        isinstance(lexpr, LiteralInt) and lexpr.value == 1
    )


def is_negative_one_lexpr(lexpr):
    """Check if an expression is negative one."""
    return (isinstance(lexpr, LiteralFloat) and lexpr.value == -1.0) or (
        isinstance(lexpr, LiteralInt) and lexpr.value == -1
    )


def float_product(factors):
    """Build product of float factors.

    Simplify ones and zeros and returning 1.0 if empty sequence.
    """
    factors = [f for f in factors if not is_one_lexpr(f)]
    if len(factors) == 0:
        return LiteralFloat(1.0)
    elif len(factors) == 1:
        return factors[0]
    else:
        for f in factors:
            if is_zero_lexpr(f):
                return f
        return Product(factors)


class DataType(Enum):
    """Representation of data types for variables in LNodes.

    These can be REAL (same type as geometry),
    SCALAR (same type as tensor), or INT (for entity indices etc.)
    """

    REAL = 0
    SCALAR = 1
    INT = 2
    BOOL = 3
    NONE = 4


def merge_dtypes(dtypes: list[DataType]):
    """Promote dtype to SCALAR or REAL if either argument matches."""
    if DataType.NONE in dtypes:
        raise ValueError(f"Invalid DataType in LNodes {dtypes}")
    if DataType.SCALAR in dtypes:
        return DataType.SCALAR
    elif DataType.REAL in dtypes:
        return DataType.REAL
    elif DataType.INT in dtypes:
        return DataType.INT
    elif DataType.BOOL in dtypes:
        return DataType.BOOL
    else:
        raise ValueError(f"Can't get dtype for operation with {dtypes}")


class LNode:
    """Base class for all AST nodes."""

    def __eq__(self, other):
        """Check for equality."""
        return NotImplemented

    def __ne__(self, other):
        """Check for inequality."""
        return NotImplemented


class LExpr(LNode):
    """Base class for all expressions.

    All subtypes should define a 'precedence' class attribute.
    """

    dtype = DataType.NONE

    def __getitem__(self, indices):
        """Get an item."""
        return ArrayAccess(self, indices)

    def __neg__(self):
        """Negate."""
        if isinstance(self, LiteralFloat):
            return LiteralFloat(-self.value)
        if isinstance(self, LiteralInt):
            return LiteralInt(-self.value)
        return Neg(self)

    def __add__(self, other):
        """Add."""
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return other
        if is_zero_lexpr(other):
            return self
        if isinstance(other, Neg):
            return Sub(self, other.arg)
        return Add(self, other)

    def __radd__(self, other):
        """Add."""
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return other
        if is_zero_lexpr(other):
            return self
        if isinstance(self, Neg):
            return Sub(other, self.arg)
        return Add(other, self)

    def __sub__(self, other):
        """Subtract."""
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return -other
        if is_zero_lexpr(other):
            return self
        if isinstance(other, Neg):
            return Add(self, other.arg)
        if isinstance(self, LiteralInt) and isinstance(other, LiteralInt):
            return LiteralInt(self.value - other.value)
        return Sub(self, other)

    def __rsub__(self, other):
        """Subtract."""
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return other
        if is_zero_lexpr(other):
            return -self
        if isinstance(self, Neg):
            return Add(other, self.arg)
        return Sub(other, self)

    def __mul__(self, other):
        """Multiply."""
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return self
        if is_zero_lexpr(other):
            return other
        if is_one_lexpr(self):
            return other
        if is_one_lexpr(other):
            return self
        if is_negative_one_lexpr(other):
            return Neg(self)
        if is_negative_one_lexpr(self):
            return Neg(other)
        if isinstance(self, LiteralInt) and isinstance(other, LiteralInt):
            return LiteralInt(self.value * other.value)
        return Mul(self, other)

    def __rmul__(self, other):
        """Multiply."""
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return self
        if is_zero_lexpr(other):
            return other
        if is_one_lexpr(self):
            return other
        if is_one_lexpr(other):
            return self
        if is_negative_one_lexpr(other):
            return Neg(self)
        if is_negative_one_lexpr(self):
            return Neg(other)
        return Mul(other, self)

    def __div__(self, other):
        """Divide."""
        other = as_lexpr(other)
        if is_zero_lexpr(other):
            raise ValueError("Division by zero!")
        if is_zero_lexpr(self):
            return self
        return Div(self, other)

    def __rdiv__(self, other):
        """Divide."""
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            raise ValueError("Division by zero!")
        if is_zero_lexpr(other):
            return other
        return Div(other, self)

    # TODO: Error check types?
    __truediv__ = __div__
    __rtruediv__ = __rdiv__
    __floordiv__ = __div__
    __rfloordiv__ = __rdiv__


class LExprOperator(LExpr):
    """Base class for all expression operators."""

    sideeffect = False


class LExprTerminal(LExpr):
    """Base class for all  expression terminals."""

    sideeffect = False


class LiteralFloat(LExprTerminal):
    """A floating point literal value."""

    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        """Initialise."""
        assert isinstance(value, (float, complex))
        self.value = value
        if isinstance(value, complex):
            self.dtype = DataType.SCALAR
        else:
            self.dtype = DataType.REAL

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, LiteralFloat) and self.value == other.value

    def __float__(self):
        """Convert to float."""
        return float(self.value)

    def __repr__(self):
        """Representation."""
        return str(self.value)


class LiteralInt(LExprTerminal):
    """An integer literal value."""

    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        """Initialise."""
        assert isinstance(value, (int, np.number))
        self.value = value
        self.dtype = DataType.INT

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, LiteralInt) and self.value == other.value

    def __hash__(self):
        """Hash."""
        return hash(self.value)

    def __repr__(self):
        """Representation."""
        return str(self.value)


class Symbol(LExprTerminal):
    """A named symbol."""

    precedence = PRECEDENCE.SYMBOL

    def __init__(self, name: str, dtype):
        """Initialise."""
        assert isinstance(name, str)
        assert name.replace("_", "").isalnum()
        self.name = name
        self.dtype = dtype

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, Symbol) and self.name == other.name

    def __hash__(self):
        """Hash."""
        return hash(self.name)

    def __repr__(self):
        """Representation."""
        return self.name


class MultiIndex(LExpr):
    """A multi-index for accessing tensors flattened in memory."""

    precedence = PRECEDENCE.SYMBOL

    def __init__(self, symbols: list, sizes: list):
        """Initialise."""
        self.dtype = DataType.INT
        self.sizes = sizes
        self.symbols = [as_lexpr(sym) for sym in symbols]
        for sym in self.symbols:
            assert sym.dtype == DataType.INT

        dim = len(sizes)
        if dim == 0:
            self.global_index: LExpr = LiteralInt(0)
        else:
            stride = [np.prod(sizes[i:]) for i in range(dim)] + [LiteralInt(1)]
            self.global_index = Sum(n * sym for n, sym in zip(stride[1:], symbols))

    @property
    def dim(self):
        """Dimension of the multi-index."""
        return len(self.sizes)

    def size(self):
        """Size of the multi-index."""
        return np.prod(self.sizes)

    def local_index(self, idx):
        """Get the local index."""
        assert idx < len(self.symbols)
        return self.symbols[idx]

    def intersection(self, other):
        """Get the intersection."""
        symbols = []
        sizes = []
        for sym, size in zip(self.symbols, self.sizes):
            if sym in other.symbols:
                i = other.symbols.index(sym)
                assert other.sizes[i] == size
                symbols.append(sym)
                sizes.append(size)
        return MultiIndex(symbols, sizes)

    def union(self, other):
        """Get the union.

        Note:
            Result may depend on order a.union(b) != b.union(a)
        """
        symbols = self.symbols.copy()
        sizes = self.sizes.copy()
        for sym, size in zip(other.symbols, other.sizes):
            if sym in symbols:
                i = symbols.index(sym)
                assert sizes[i] == size
            else:
                symbols.append(sym)
                sizes.append(size)
        return MultiIndex(symbols, sizes)

    def difference(self, other):
        """Get the difference."""
        symbols = []
        sizes = []
        for idx, size in zip(self.symbols, self.sizes):
            if idx not in other.symbols:
                symbols.append(idx)
                sizes.append(size)
        return MultiIndex(symbols, sizes)

    def __hash__(self):
        """Hash."""
        return hash(self.global_index.__repr__)


class PrefixUnaryOp(LExprOperator):
    """Base class for unary operators."""

    def __init__(self, arg):
        """Initialise."""
        self.arg = as_lexpr(arg)

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, type(self)) and self.arg == other.arg


class BinOp(LExprOperator):
    """A binary operator."""

    def __init__(self, lhs, rhs):
        """Initialise."""
        self.lhs = as_lexpr(lhs)
        self.rhs = as_lexpr(rhs)

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, type(self)) and self.lhs == other.lhs and self.rhs == other.rhs

    def __hash__(self):
        """Hash."""
        return hash(self.lhs) + hash(self.rhs)

    def __repr__(self):
        """Representation."""
        return f"({self.lhs} {self.op} {self.rhs})"


class ArithmeticBinOp(BinOp):
    """An artithmetic binary operator."""

    def __init__(self, lhs, rhs):
        """Initialise."""
        self.lhs = as_lexpr(lhs)
        self.rhs = as_lexpr(rhs)
        self.dtype = merge_dtypes([self.lhs.dtype, self.rhs.dtype])


class NaryOp(LExprOperator):
    """Base class for special n-ary operators."""

    op = ""

    def __init__(self, args):
        """Initialise."""
        self.args = [as_lexpr(arg) for arg in args]
        self.dtype = self.args[0].dtype
        for arg in self.args:
            self.dtype = merge_dtypes([self.dtype, arg.dtype])

    def __eq__(self, other):
        """Check equality."""
        return (
            isinstance(other, type(self))
            and len(self.args) == len(other.args)
            and all(a == b for a, b in zip(self.args, other.args))
        )

    def __repr__(self) -> str:
        """Representation."""
        return f"{self.op} ".join(f"{i} " for i in self.args)

    def __hash__(self):
        """Hash."""
        return hash(tuple(self.args))


class Neg(PrefixUnaryOp):
    """Negation operator."""

    precedence = PRECEDENCE.NEG
    op = "-"

    def __init__(self, arg):
        """Initialise."""
        self.arg = as_lexpr(arg)
        self.dtype = self.arg.dtype


class Not(PrefixUnaryOp):
    """Not operator."""

    precedence = PRECEDENCE.NOT
    op = "!"


class Add(ArithmeticBinOp):
    """Add operator."""

    precedence = PRECEDENCE.ADD
    op = "+"


class Sub(ArithmeticBinOp):
    """Subtract operator."""

    precedence = PRECEDENCE.SUB
    op = "-"


class Mul(ArithmeticBinOp):
    """Multiply operator."""

    precedence = PRECEDENCE.MUL
    op = "*"


class Div(ArithmeticBinOp):
    """Division operator."""

    precedence = PRECEDENCE.DIV
    op = "/"


class EQ(BinOp):
    """Equality operator."""

    precedence = PRECEDENCE.EQ
    op = "=="


class NE(BinOp):
    """Inequality operator."""

    precedence = PRECEDENCE.NE
    op = "!="


class LT(BinOp):
    """Less than operator."""

    precedence = PRECEDENCE.LT
    op = "<"


class GT(BinOp):
    """Greater than operator."""

    precedence = PRECEDENCE.GT
    op = ">"


class LE(BinOp):
    """Less than or equal to operator."""

    precedence = PRECEDENCE.LE
    op = "<="


class GE(BinOp):
    """Greater than or equal to operator."""

    precedence = PRECEDENCE.GE
    op = ">="


class And(BinOp):
    """And operator."""

    precedence = PRECEDENCE.AND
    op = "&&"


class Or(BinOp):
    """Or operator."""

    precedence = PRECEDENCE.OR
    op = "||"


class Sum(NaryOp):
    """Sum of any number of operands."""

    precedence = PRECEDENCE.ADD
    op = "+"


class Product(NaryOp):
    """Product of any number of operands."""

    precedence = PRECEDENCE.MUL
    op = "*"


class MathFunction(LExprOperator):
    """A Math Function, with any arguments."""

    precedence = PRECEDENCE.HIGHEST

    def __init__(self, func, args):
        """Initialise."""
        self.function = func
        self.args = [as_lexpr(arg) for arg in args]
        self.dtype = self.args[0].dtype

    def __eq__(self, other):
        """Check equality."""
        return (
            isinstance(other, type(self))
            and self.function == other.function
            and len(self.args) == len(other.args)
            and all(a == b for a, b in zip(self.args, other.args))
        )


class AssignOp(BinOp):
    """Base class for assignment operators."""

    precedence = PRECEDENCE.ASSIGN
    sideeffect = True

    def __init__(self, lhs, rhs):
        """Initialise."""
        assert isinstance(lhs, LNode)
        BinOp.__init__(self, lhs, rhs)


class Assign(AssignOp):
    """Assign operator."""

    op = "="


class AssignAdd(AssignOp):
    """Assign add operator."""

    op = "+="


class AssignSub(AssignOp):
    """Assign subtract operator."""

    op = "-="


class AssignMul(AssignOp):
    """Assign multiply operator."""

    op = "*="


class AssignDiv(AssignOp):
    """Assign division operator."""

    op = "/="


class ArrayAccess(LExprOperator):
    """Array access."""

    precedence = PRECEDENCE.SUBSCRIPT

    def __init__(self, array, indices):
        """Initialise."""
        # Typecheck array argument
        if isinstance(array, Symbol):
            self.array = array
            self.dtype = array.dtype
        elif isinstance(array, ArrayDecl):
            self.array = array.symbol
            self.dtype = array.symbol.dtype
        else:
            raise ValueError(f"Unexpected array type {type(array).__name__}")

        # Allow expressions or literals as indices
        if not isinstance(indices, (list, tuple)):
            indices = (indices,)
        self.indices = tuple(as_lexpr(i) for i in indices)

        # Early error checking for negative array dimensions
        if any(isinstance(i, int) and i < 0 for i in self.indices):
            raise ValueError("Index value < 0.")

        # Additional dimension checks possible if we get an ArrayDecl instead of just a name
        if isinstance(array, ArrayDecl):
            if len(self.indices) != len(array.sizes):
                raise ValueError("Invalid number of indices.")
            ints = (int, LiteralInt)
            if any(
                (isinstance(i, ints) and isinstance(d, ints) and int(i) >= int(d))
                for i, d in zip(self.indices, array.sizes)
            ):
                raise ValueError("Index value >= array dimension.")

    def __getitem__(self, indices):
        """Handle nested expr[i][j]."""
        if isinstance(indices, list):
            indices = tuple(indices)
        elif not isinstance(indices, tuple):
            indices = (indices,)
        return ArrayAccess(self.array, self.indices + indices)

    def __eq__(self, other):
        """Check equality."""
        return (
            isinstance(other, type(self))
            and self.array == other.array
            and self.indices == other.indices
        )

    def __hash__(self):
        """Hash."""
        return hash(self.array)

    def __repr__(self):
        """Representation."""
        return str(self.array) + "[" + ", ".join(str(i) for i in self.indices) + "]"


class Conditional(LExprOperator):
    """Conditional."""

    precedence = PRECEDENCE.CONDITIONAL

    def __init__(self, condition, true, false):
        """Initialise."""
        self.condition = as_lexpr(condition)
        self.true = as_lexpr(true)
        self.false = as_lexpr(false)
        self.dtype = merge_dtypes([self.true.dtype, self.false.dtype])

    def __eq__(self, other):
        """Check equality."""
        return (
            isinstance(other, type(self))
            and self.condition == other.condition
            and self.true == other.true
            and self.false == other.false
        )


def as_lexpr(node):
    """Typechecks and wraps an object as a valid LExpr.

    Accepts LExpr nodes, treats int and float as literals.

    """
    if isinstance(node, LExpr):
        return node
    elif isinstance(node, numbers.Integral):
        return LiteralInt(node)
    elif isinstance(node, numbers.Real):
        return LiteralFloat(node)
    else:
        raise RuntimeError(f"Unexpected LExpr type {type(node)}:\n{node}")


class Statement(LNode):
    """Make an expression into a statement."""

    def __init__(self, expr):
        """Initialise."""
        self.expr = as_lexpr(expr)

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, type(self)) and self.expr == other.expr

    def __hash__(self) -> int:
        """Hash."""
        return hash(self.expr)


def as_statement(node):
    """Perform type checking on node and wrap in a suitable statement type if necessary."""
    if isinstance(node, StatementList) and len(node.statements) == 1:
        # Cleans up the expression tree a bit
        return node.statements[0]
    elif isinstance(node, Statement):
        # No-op
        return node
    elif isinstance(node, LExprOperator):
        if node.sideeffect:
            # Special case for using assignment expressions as statements
            return Statement(node)
        else:
            raise RuntimeError(
                f"Trying to create a statement of lexprOperator type {type(node)}:\n{node}"
            )

    elif isinstance(node, list):
        # Convenience case for list of statements
        if len(node) == 1:
            # Cleans up the expression tree a bit
            return as_statement(node[0])
        else:
            return StatementList(node)
    elif isinstance(node, Section):
        return node
    else:
        raise RuntimeError(f"Unexpected Statement type {type(node)}:\n{node}")


class Annotation(Enum):
    """Annotation."""

    fuse = 1  # fuse loops in section
    unroll = 2  # unroll loop in section
    licm = 3  # loop invariant code motion
    factorize = 4  # apply sum factorization


class Declaration(Statement):
    """Base class for all declarations."""

    def __init__(self, symbol):
        """Initialise."""
        self.symbol = symbol

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, type(self)) and self.symbol == other.symbol


def is_declaration(node) -> bool:
    """Check if a node is a declaration."""
    return isinstance(node, VariableDecl) or isinstance(node, ArrayDecl)


class Section(LNode):
    """A section of code with a name and a list of statements."""

    def __init__(
        self,
        name: str,
        statements: list[LNode],
        declarations: Sequence[Declaration],
        input: Optional[list[Symbol]] = None,
        output: Optional[list[Symbol]] = None,
        annotations: Optional[list[Annotation]] = None,
    ):
        """Initialise."""
        self.name = name
        self.statements = [as_statement(st) for st in statements]
        self.annotations = annotations or []
        self.input = input or []
        self.declarations = declarations or []
        self.output = output or []

        for decl in self.declarations:
            assert is_declaration(decl)
            if decl.symbol not in self.output:
                self.output.append(decl.symbol)

    def __eq__(self, other):
        """Check equality."""
        attributes = ("name", "input", "output", "annotations", "statements")
        return isinstance(other, type(self)) and all(
            getattr(self, name) == getattr(self, name) for name in attributes
        )


class StatementList(LNode):
    """A simple sequence of statements."""

    def __init__(self, statements):
        """Initialise."""
        self.statements = [as_statement(st) for st in statements]

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, type(self)) and self.statements == other.statements

    def __hash__(self) -> int:
        """Hash."""
        return hash(tuple(self.statements))

    def __repr__(self):
        """Representation."""
        return f"StatementList({self.statements})"


class Comment(Statement):
    """Line comment(s) used for annotating the generated code with human readable remarks."""

    def __init__(self, comment):
        """Initialise."""
        assert isinstance(comment, str)
        self.comment = comment

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, type(self)) and self.comment == other.comment


def commented_code_list(code, comments):
    """Add comment to code list if the list is not empty."""
    if isinstance(code, LNode):
        code = [code]
    assert isinstance(code, list)
    if code:
        if not isinstance(comments, (list, tuple)):
            comments = [comments]
        comments = [Comment(c) for c in comments]
        code = comments + code
    return code


# Type and variable declarations


class VariableDecl(Declaration):
    """Declare a variable, optionally define initial value."""

    def __init__(self, symbol, value=None):
        """Initialise."""
        assert isinstance(symbol, Symbol)
        assert symbol.dtype is not None
        self.symbol = symbol

        if value is not None:
            value = as_lexpr(value)
        self.value = value

    def __eq__(self, other):
        """Check equality."""
        return (
            isinstance(other, type(self))
            and self.typename == other.typename
            and self.symbol == other.symbol
            and self.value == other.value
        )


class ArrayDecl(Declaration):
    """A declaration or definition of an array.

    Note that just setting values=0 is sufficient to initialize the
    entire array to zero.

    Otherwise use nested lists of lists to represent multidimensional
    array values to initialize to.

    """

    def __init__(self, symbol, sizes=None, values=None, const=False):
        """Initialise."""
        assert isinstance(symbol, Symbol)
        self.symbol = symbol
        assert symbol.dtype

        if sizes is None:
            assert values is not None
            sizes = values.shape
        if isinstance(sizes, int):
            sizes = (sizes,)
        self.sizes = tuple(sizes)

        if values is None:
            assert sizes is not None

        # NB! No type checking, assuming nested lists of literal values. Not applying as_lexpr.
        if isinstance(values, (list, tuple)):
            self.values = np.asarray(values)
        else:
            self.values = values

        self.const = const
        self.dtype = symbol.dtype

    def __eq__(self, other):
        """Check equality."""
        attributes = ("dtype", "symbol", "sizes", "values")
        return isinstance(other, type(self)) and all(
            getattr(self, name) == getattr(self, name) for name in attributes
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(self.symbol)


def is_simple_inner_loop(code):
    """Check if code is a simple inner loop."""
    if isinstance(code, ForRange) and is_simple_inner_loop(code.body):
        return True
    if isinstance(code, Statement) and isinstance(code.expr, AssignOp):
        return True
    return False


def depth(code) -> int:
    """Get depth of code."""
    if isinstance(code, ForRange):
        return 1 + depth(code.body)
    if isinstance(code, StatementList):
        return max([depth(c) for c in code.statements])
    return 0


class ForRange(Statement):
    """Slightly higher-level for loop assuming incrementing an index over a range."""

    def __init__(self, index, begin, end, body):
        """Initialise."""
        assert isinstance(index, Symbol) or isinstance(index, MultiIndex)
        self.index = index
        self.begin = as_lexpr(begin)
        self.end = as_lexpr(end)
        assert isinstance(body, list)
        self.body = StatementList(body)

    def as_tuple(self):
        """Convert to a tuple."""
        return (self.index, self.begin, self.end, self.body)

    def __eq__(self, other):
        """Check equality."""
        attributes = ("index", "begin", "end", "body")
        return isinstance(other, type(self)) and all(
            getattr(self, name) == getattr(self, name) for name in attributes
        )

    def __hash__(self) -> int:
        """Hash."""
        return hash(self.as_tuple())


def _math_function(op, *args):
    """Get a math function."""
    name = op._ufl_handler_name_
    dtype = args[0].dtype
    if name in ("conj", "real") and dtype == DataType.REAL:
        assert len(args) == 1
        return args[0]
    if name == "imag" and dtype == DataType.REAL:
        assert len(args) == 1
        return LiteralFloat(0.0)
    return MathFunction(name, args)


# Lookup table for handler to call when the ufl_to_lnodes method (below) is
# called, depending on the first argument type.
_ufl_call_lookup = {
    ufl.constantvalue.IntValue: lambda x: LiteralInt(int(x)),
    ufl.constantvalue.FloatValue: lambda x: LiteralFloat(float(x)),
    ufl.constantvalue.ComplexValue: lambda x: LiteralFloat(x.value()),
    ufl.constantvalue.Zero: lambda x: LiteralFloat(0.0),
    ufl.algebra.Product: lambda x, a, b: a * b,
    ufl.algebra.Sum: lambda x, a, b: a + b,
    ufl.algebra.Division: lambda x, a, b: a / b,
    ufl.algebra.Abs: _math_function,
    ufl.algebra.Power: _math_function,
    ufl.algebra.Real: _math_function,
    ufl.algebra.Imag: _math_function,
    ufl.algebra.Conj: _math_function,
    ufl.classes.GT: lambda x, a, b: GT(a, b),
    ufl.classes.GE: lambda x, a, b: GE(a, b),
    ufl.classes.EQ: lambda x, a, b: EQ(a, b),
    ufl.classes.NE: lambda x, a, b: NE(a, b),
    ufl.classes.LT: lambda x, a, b: LT(a, b),
    ufl.classes.LE: lambda x, a, b: LE(a, b),
    ufl.classes.AndCondition: lambda x, a, b: And(a, b),
    ufl.classes.OrCondition: lambda x, a, b: Or(a, b),
    ufl.classes.NotCondition: lambda x, a: Not(a),
    ufl.classes.Conditional: lambda x, c, t, f: Conditional(c, t, f),
    ufl.classes.MinValue: _math_function,
    ufl.classes.MaxValue: _math_function,
    ufl.mathfunctions.Sqrt: _math_function,
    ufl.mathfunctions.Ln: _math_function,
    ufl.mathfunctions.Exp: _math_function,
    ufl.mathfunctions.Cos: _math_function,
    ufl.mathfunctions.Sin: _math_function,
    ufl.mathfunctions.Tan: _math_function,
    ufl.mathfunctions.Cosh: _math_function,
    ufl.mathfunctions.Sinh: _math_function,
    ufl.mathfunctions.Tanh: _math_function,
    ufl.mathfunctions.Acos: _math_function,
    ufl.mathfunctions.Asin: _math_function,
    ufl.mathfunctions.Atan: _math_function,
    ufl.mathfunctions.Erf: _math_function,
    ufl.mathfunctions.Atan2: _math_function,
    ufl.mathfunctions.MathFunction: _math_function,
    ufl.mathfunctions.BesselJ: _math_function,
    ufl.mathfunctions.BesselY: _math_function,
}


def ufl_to_lnodes(operator, *args):
    """Call appropriate handler, depending on the type of operator."""
    optype = type(operator)
    if optype in _ufl_call_lookup:
        return _ufl_call_lookup[optype](operator, *args)
    else:
        raise RuntimeError(f"Missing lookup for expr type {optype}.")


def create_nested_for_loops(indices: list[MultiIndex], body):
    """Create nested for loops over list of indices.

    The depth of the nested for loops is equal to the sub-indices for all
    MultiIndex combined.
    """
    ranges = [r for idx in indices for r in idx.sizes]
    indices = [idx.local_index(i) for idx in indices for i in range(len(idx.sizes))]
    depth = len(ranges)
    for i in reversed(range(depth)):
        body = ForRange(indices[i], 0, ranges[i], body=[body])
    return body
