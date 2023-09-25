# Copyright (C) 2013-2023 Martin Sandve Alnæs, Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numbers
import ufl
import numpy as np
from enum import Enum


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


"""LNodes is intended as a minimal generic language description.
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


def is_zero_lexpr(lexpr):
    return (isinstance(lexpr, LiteralFloat) and lexpr.value == 0.0) or (
        isinstance(lexpr, LiteralInt) and lexpr.value == 0
    )


def is_one_lexpr(lexpr):
    return (isinstance(lexpr, LiteralFloat) and lexpr.value == 1.0) or (
        isinstance(lexpr, LiteralInt) and lexpr.value == 1
    )


def is_negative_one_lexpr(lexpr):
    return (isinstance(lexpr, LiteralFloat) and lexpr.value == -1.0) or (
        isinstance(lexpr, LiteralInt) and lexpr.value == -1
    )


def float_product(factors):
    """Build product of float factors, simplifying ones and zeros and returning 1.0 if empty sequence."""
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
    NONE = 3


def merge_dtypes(dtype0, dtype1):
    # Promote dtype to SCALAR or REAL if either argument matches
    if DataType.NONE in (dtype0, dtype1):
        raise ValueError(f"Invalid DataType in LNodes {dtype0, dtype1}")
    if DataType.SCALAR in (dtype0, dtype1):
        return DataType.SCALAR
    elif DataType.REAL in (dtype0, dtype1):
        return DataType.REAL
    elif (dtype0 == DataType.INT and dtype1 == DataType.INT):
        return DataType.INT
    else:
        raise ValueError(f"Can't get dtype for binary operation with {dtype0, dtype1}")


class LNode(object):
    """Base class for all AST nodes."""

    def __eq__(self, other):
        name = self.__class__.__name__
        raise NotImplementedError("Missing implementation of __eq__ in " + name)

    def __ne__(self, other):
        return not self.__eq__(other)


class LExpr(LNode):
    """Base class for all expressions.

    All subtypes should define a 'precedence' class attribute.
    """

    dtype = DataType.NONE

    def __getitem__(self, indices):
        return ArrayAccess(self, indices)

    def __neg__(self):
        if isinstance(self, LiteralFloat):
            return LiteralFloat(-self.value)
        if isinstance(self, LiteralInt):
            return LiteralInt(-self.value)
        return Neg(self)

    def __add__(self, other):
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return other
        if is_zero_lexpr(other):
            return self
        if isinstance(other, Neg):
            return Sub(self, other.arg)
        return Add(self, other)

    def __radd__(self, other):
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return other
        if is_zero_lexpr(other):
            return self
        if isinstance(self, Neg):
            return Sub(other, self.arg)
        return Add(other, self)

    def __sub__(self, other):
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
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            return other
        if is_zero_lexpr(other):
            return -self
        if isinstance(self, Neg):
            return Add(other, self.arg)
        return Sub(other, self)

    def __mul__(self, other):
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
        other = as_lexpr(other)
        if is_zero_lexpr(other):
            raise ValueError("Division by zero!")
        if is_zero_lexpr(self):
            return self
        return Div(self, other)

    def __rdiv__(self, other):
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


# LExprTerminal types


class LiteralFloat(LExprTerminal):
    """A floating point literal value."""

    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (float, complex))
        self.value = value
        if isinstance(value, complex):
            self.dtype = DataType.SCALAR
        else:
            self.dtype = DataType.REAL

    def __eq__(self, other):
        return isinstance(other, LiteralFloat) and self.value == other.value

    def __float__(self):
        return float(self.value)

    def __repr__(self):
        return str(self.value)


class LiteralInt(LExprTerminal):
    """An integer literal value."""

    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (int, np.number))
        self.value = value
        self.dtype = DataType.INT

    def __eq__(self, other):
        return isinstance(other, LiteralInt) and self.value == other.value

    def __hash__(self):
        return hash(self.value)

    def __repr__(self):
        return str(self.value)


class Symbol(LExprTerminal):
    """A named symbol."""

    precedence = PRECEDENCE.SYMBOL

    def __init__(self, name: str, dtype):
        assert isinstance(name, str)
        assert name.replace("_", "").isalnum()
        self.name = name
        self.dtype = dtype

    def __eq__(self, other):
        return isinstance(other, Symbol) and self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return self.name


class MultiIndex(LExpr):
    """A multi-index for accessing tensors flattened in memory."""

    def __init__(self, symbols: list, sizes: list):
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

    def size(self):
        return np.prod(self.sizes)

    def local_index(self, idx):
        assert idx < len(self.symbols)
        return self.symbols[idx]

    def intersection(self, other):
        symbols = []
        sizes = []
        for (sym, size) in zip(self.symbols, self.sizes):
            if sym in other.symbols:
                i = other.symbols.index(sym)
                assert other.sizes[i] == size
                symbols.append(sym)
                sizes.append(size)
        return MultiIndex(symbols, sizes)

    def union(self, other):
        # NB result may depend on order a.union(b) != b.union(a)
        symbols = self.symbols.copy()
        sizes = self.sizes.copy()
        for (sym, size) in zip(other.symbols, other.sizes):
            if sym in symbols:
                i = symbols.index(sym)
                assert sizes[i] == size
            else:
                symbols.append(sym)
                sizes.append(size)
        return MultiIndex(symbols, sizes)

    def difference(self, other):
        symbols = []
        sizes = []
        for (idx, size) in zip(self.symbols, self.sizes):
            if idx not in other.symbols:
                symbols.append(idx)
                sizes.append(size)
        return MultiIndex(symbols, sizes)

    def __hash__(self):
        return hash(self.global_idx)


class PrefixUnaryOp(LExprOperator):
    """Base class for unary operators."""

    def __init__(self, arg):
        self.arg = as_lexpr(arg)

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.arg == other.arg


class BinOp(LExprOperator):
    def __init__(self, lhs, rhs):
        self.lhs = as_lexpr(lhs)
        self.rhs = as_lexpr(rhs)

    def __eq__(self, other):
        return (
            isinstance(other, type(self))
            and self.lhs == other.lhs
            and self.rhs == other.rhs
        )

    def __hash__(self):
        return hash(self.lhs) + hash(self.rhs)

    def __repr__(self):
        return f"({self.lhs} {self.op} {self.rhs})"


class ArithmeticBinOp(BinOp):
    def __init__(self, lhs, rhs):
        self.lhs = as_lexpr(lhs)
        self.rhs = as_lexpr(rhs)
        self.dtype = merge_dtypes(self.lhs.dtype, self.rhs.dtype)


class NaryOp(LExprOperator):
    """Base class for special n-ary operators."""

    op = ""

    def __init__(self, args):
        self.args = [as_lexpr(arg) for arg in args]

    def __eq__(self, other):
        return (
            isinstance(other, type(self))
            and len(self.args) == len(other.args)
            and all(a == b for a, b in zip(self.args, other.args))
        )

    def __repr__(self) -> str:
        return f"{self.op} ".join(f"{i} " for i in self.args)


class Neg(PrefixUnaryOp):
    precedence = PRECEDENCE.NEG
    op = "-"

    def __init__(self, arg):
        self.arg = as_lexpr(arg)
        self.dtype = self.arg.dtype


class Not(PrefixUnaryOp):
    precedence = PRECEDENCE.NOT
    op = "!"


# Binary operators
# Arithmetic operators preserve the dtype of their operands
# The other operations (logical) do not need a dtype

class Add(ArithmeticBinOp):
    precedence = PRECEDENCE.ADD
    op = "+"


class Sub(ArithmeticBinOp):
    precedence = PRECEDENCE.SUB
    op = "-"


class Mul(ArithmeticBinOp):
    precedence = PRECEDENCE.MUL
    op = "*"


class Div(ArithmeticBinOp):
    precedence = PRECEDENCE.DIV
    op = "/"


class EQ(BinOp):
    precedence = PRECEDENCE.EQ
    op = "=="


class NE(BinOp):
    precedence = PRECEDENCE.NE
    op = "!="


class LT(BinOp):
    precedence = PRECEDENCE.LT
    op = "<"


class GT(BinOp):
    precedence = PRECEDENCE.GT
    op = ">"


class LE(BinOp):
    precedence = PRECEDENCE.LE
    op = "<="


class GE(BinOp):
    precedence = PRECEDENCE.GE
    op = ">="


class And(BinOp):
    precedence = PRECEDENCE.AND
    op = "&&"


class Or(BinOp):
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
        self.function = func
        self.args = [as_lexpr(arg) for arg in args]
        self.dtype = self.args[0].dtype

    def __eq__(self, other):
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
        assert isinstance(lhs, LNode)
        BinOp.__init__(self, lhs, rhs)


class Assign(AssignOp):
    op = "="


class AssignAdd(AssignOp):
    op = "+="


class AssignSub(AssignOp):
    op = "-="


class AssignMul(AssignOp):
    op = "*="


class AssignDiv(AssignOp):
    op = "/="


class ArrayAccess(LExprOperator):
    precedence = PRECEDENCE.SUBSCRIPT

    def __init__(self, array, indices):
        # Typecheck array argument
        if isinstance(array, Symbol):
            self.array = array
            self.dtype = array.dtype
        elif isinstance(array, ArrayDecl):
            self.array = array.symbol
            self.dtype = array.symbol.dtype
        else:
            raise ValueError("Unexpected array type %s." % (type(array).__name__,))

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
        return (
            isinstance(other, type(self))
            and self.array == other.array
            and self.indices == other.indices
        )

    def __hash__(self):
        return hash(self.array)

    def __repr__(self):
        return str(self.array) + "[" + ", ".join(str(i) for i in self.indices) + "]"


class Conditional(LExprOperator):
    precedence = PRECEDENCE.CONDITIONAL

    def __init__(self, condition, true, false):
        self.condition = as_lexpr(condition)
        self.true = as_lexpr(true)
        self.false = as_lexpr(false)
        self.dtype = merge_dtypes(self.true.dtype, self.false.dtype)

    def __eq__(self, other):
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
        raise RuntimeError("Unexpected LExpr type %s:\n%s" % (type(node), str(node)))


class Statement(LNode):
    """Make an expression into a statement."""

    is_scoped = False

    def __init__(self, expr):
        self.expr = as_lexpr(expr)

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.expr == other.expr


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
                "Trying to create a statement of lexprOperator type %s:\n%s"
                % (type(node), str(node))
            )
    elif isinstance(node, list):
        # Convenience case for list of statements
        if len(node) == 1:
            # Cleans up the expression tree a bit
            return as_statement(node[0])
        else:
            return StatementList(node)
    else:
        raise RuntimeError(
            "Unexpected Statement type %s:\n%s" % (type(node), str(node))
        )


class StatementList(LNode):
    """A simple sequence of statements. No new scopes are introduced."""

    def __init__(self, statements):
        self.statements = [as_statement(st) for st in statements]

    @property
    def is_scoped(self):
        return all(st.is_scoped for st in self.statements)

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.statements == other.statements


class Comment(Statement):
    """Line comment(s) used for annotating the generated code with human readable remarks."""

    is_scoped = True

    def __init__(self, comment):
        assert isinstance(comment, str)
        self.comment = comment

    def __eq__(self, other):
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


class VariableDecl(Statement):
    """Declare a variable, optionally define initial value."""

    is_scoped = False

    def __init__(self, symbol, value=None):

        assert isinstance(symbol, Symbol)
        assert symbol.dtype is not None
        self.symbol = symbol

        if value is not None:
            value = as_lexpr(value)
        self.value = value

    def __eq__(self, other):
        return (
            isinstance(other, type(self))
            and self.typename == other.typename
            and self.symbol == other.symbol
            and self.value == other.value
        )


class ArrayDecl(Statement):
    """A declaration or definition of an array.

    Note that just setting values=0 is sufficient to initialize the
    entire array to zero.

    Otherwise use nested lists of lists to represent multidimensional
    array values to initialize to.

    """

    is_scoped = False

    def __init__(self, symbol, sizes=None, values=None, const=False):
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

    def __eq__(self, other):
        attributes = ("typename", "symbol", "sizes", "values")
        return isinstance(other, type(self)) and all(
            getattr(self, name) == getattr(self, name) for name in attributes
        )


def is_simple_inner_loop(code):
    if isinstance(code, ForRange) and is_simple_inner_loop(code.body):
        return True
    if isinstance(code, Statement) and isinstance(code.expr, AssignOp):
        return True
    return False


class ForRange(Statement):
    """Slightly higher-level for loop assuming incrementing an index over a range."""

    is_scoped = True

    def __init__(self, index, begin, end, body):
        assert isinstance(index, Symbol)
        self.index = index
        self.begin = as_lexpr(begin)
        self.end = as_lexpr(end)
        assert isinstance(body, list)
        self.body = StatementList(body)

    def __eq__(self, other):
        attributes = ("index", "begin", "end", "body")
        return isinstance(other, type(self)) and all(
            getattr(self, name) == getattr(self, name) for name in attributes
        )


def _math_function(op, *args):
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
    ufl.mathfunctions.BesselY: _math_function}


def ufl_to_lnodes(operator, *args):
    # Call appropriate handler, depending on the type of operator
    optype = type(operator)
    if optype in _ufl_call_lookup:
        return _ufl_call_lookup[optype](operator, *args)
    else:
        raise RuntimeError(f"Missing lookup for expr type {optype}.")
