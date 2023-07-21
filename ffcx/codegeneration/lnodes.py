# Copyright (C) 2013-2023 Martin Sandve Alnæs, Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numbers

import numpy as np


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
    MOD = 4

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


"""Intended as a minimal generic language description. 
Formatting is done later, depending on the target language.

Supported:
 Floating point (and complex) and integer variables and multidimensional arrays
 Range loops
 Simple arithmetic, +-*/
 Math operations
 Comments
Not supported:
 Pointers
 Function Calls
 Flow control (if, switch, while)
 Booleans
 Logical operations
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

    # TODO: Error check types? Can't do that exactly as symbols here have no type.
    __truediv__ = __div__
    __rtruediv__ = __rdiv__
    __floordiv__ = __div__
    __rfloordiv__ = __rdiv__

    def __mod__(self, other):
        other = as_lexpr(other)
        if is_zero_lexpr(other):
            raise ValueError("Division by zero!")
        if is_zero_lexpr(self):
            return self
        return Mod(self, other)

    def __rmod__(self, other):
        other = as_lexpr(other)
        if is_zero_lexpr(self):
            raise ValueError("Division by zero!")
        if is_zero_lexpr(other):
            return other
        return Mod(other, self)


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
        assert isinstance(value, (float, complex, int, np.number))
        self.value = value

    def __eq__(self, other):
        return isinstance(other, LiteralFloat) and self.value == other.value

    def __float__(self):
        return float(self.value)


class LiteralInt(LExprTerminal):
    """An integer literal value."""

    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (int, np.number))
        self.value = value

    def __eq__(self, other):
        return isinstance(other, LiteralInt) and self.value == other.value

    def __hash__(self):
        return hash(self.value)


class Symbol(LExprTerminal):
    """A named symbol."""

    precedence = PRECEDENCE.SYMBOL

    def __init__(self, name):
        assert isinstance(name, str)
        self.name = name

    def __eq__(self, other):
        return isinstance(other, Symbol) and self.name == other.name

    def __hash__(self):
        return hash(self.name)


class UnaryOp(LExprOperator):
    """Base class for unary operators."""

    def __init__(self, arg):
        self.arg = as_lexpr(arg)

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.arg == other.arg


class PrefixUnaryOp(UnaryOp):
    """Base class for prefix unary operators."""

    def __eq__(self, other):
        return isinstance(other, type(self))


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


class NaryOp(LExprOperator):
    """Base class for special n-ary operators."""

    def __init__(self, args):
        self.args = [as_lexpr(arg) for arg in args]

    def __eq__(self, other):
        return (
            isinstance(other, type(self))
            and len(self.args) == len(other.args)
            and all(a == b for a, b in zip(self.args, other.args))
        )


class Neg(PrefixUnaryOp):
    precedence = PRECEDENCE.NEG
    op = "-"


class Not(PrefixUnaryOp):
    precedence = PRECEDENCE.NOT
    op = "!"


# lexpr binary operators


class Add(BinOp):
    precedence = PRECEDENCE.ADD
    op = "+"


class Sub(BinOp):
    precedence = PRECEDENCE.SUB
    op = "-"


class Mul(BinOp):
    precedence = PRECEDENCE.MUL
    op = "*"


class Div(BinOp):
    precedence = PRECEDENCE.DIV
    op = "/"


class Mod(BinOp):
    precedence = PRECEDENCE.MOD
    op = "%"


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
    """A Math Function, with any arguments"""

    precedence = PRECEDENCE.HIGHEST

    def __init__(self, func, args):
        self.function = func
        self.args = [as_lexpr(arg) for arg in args]

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
        BinOp.__init__(self, as_lexpr_or_string_symbol(lhs), rhs)


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


class FlattenedArray(object):
    """Syntax carrying object only, will get translated on __getitem__ to ArrayAccess."""

    def __init__(self, array, dummy=None, dims=None, strides=None, offset=None):
        assert dummy is None, "Please use keyword arguments for strides or dims."

        # Typecheck array argument
        if isinstance(array, ArrayDecl):
            self.array = array.symbol
        elif isinstance(array, Symbol):
            self.array = array
        else:
            assert isinstance(array, str)
            self.array = Symbol(array)

        # Allow expressions or literals as strides or dims and offset
        if strides is None:
            assert dims is not None, "Please provide either strides or dims."
            assert isinstance(dims, (list, tuple))
            dims = tuple(as_lexpr(i) for i in dims)
            self.dims = dims
            n = len(dims)
            literal_one = LiteralInt(1)
            strides = [literal_one] * n
            for i in range(n - 2, -1, -1):
                s = strides[i + 1]
                d = dims[i + 1]
                if d == literal_one:
                    strides[i] = s
                elif s == literal_one:
                    strides[i] = d
                else:
                    strides[i] = d * s
        else:
            self.dims = None
            assert isinstance(strides, (list, tuple))
            strides = tuple(as_lexpr(i) for i in strides)
        self.strides = strides
        self.offset = None if offset is None else as_lexpr(offset)

    def __getitem__(self, indices):
        if not isinstance(indices, (list, tuple)):
            indices = (indices,)
        n = len(indices)
        if n == 0:
            # Handle scalar case, allowing dims=() and indices=() for A[0]
            if len(self.strides) != 0:
                raise ValueError("Empty indices for nonscalar array.")
            flat = LiteralInt(0)
        else:
            i, s = (indices[0], self.strides[0])
            literal_one = LiteralInt(1)
            flat = i if s == literal_one else s * i
            if self.offset is not None:
                flat = self.offset + flat
            for i, s in zip(indices[1:n], self.strides[1:n]):
                flat = flat + (i if s == literal_one else s * i)
        # Delay applying ArrayAccess until we have all indices
        if n == len(self.strides):
            return ArrayAccess(self.array, flat)
        else:
            return FlattenedArray(self.array, strides=self.strides[n:], offset=flat)


class ArrayAccess(LExprOperator):
    precedence = PRECEDENCE.SUBSCRIPT

    def __init__(self, array, indices):
        # Typecheck array argument
        if isinstance(array, Symbol):
            self.array = array
        elif isinstance(array, ArrayDecl):
            self.array = array.symbol
        else:
            raise ValueError("Unexpected array type %s." % (type(array).__name__,))

        # Allow expressions or literals as indices
        if not isinstance(indices, (list, tuple)):
            indices = (indices,)
        self.indices = tuple(as_lexpr_or_string_symbol(i) for i in indices)

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


class Conditional(LExprOperator):
    precedence = PRECEDENCE.CONDITIONAL

    def __init__(self, condition, true, false):
        self.condition = as_lexpr(condition)
        self.true = as_lexpr(true)
        self.false = as_lexpr(false)

    def __eq__(self, other):
        return (
            isinstance(other, type(self))
            and self.condition == other.condition
            and self.true == other.true
            and self.false == other.false
        )


def _is_zero_valued(values):
    if isinstance(values, (numbers.Integral, LiteralInt)):
        return int(values) == 0
    elif isinstance(values, (numbers.Number, LiteralFloat)):
        return float(values) == 0.0
    else:
        return np.count_nonzero(values) == 0


def as_lexpr(node):
    """Typechecks and wraps an object as a valid LExpr.

    Accepts LExpr nodes, treats int and float as literals, and treats a
    string as a symbol.

    """
    if isinstance(node, LExpr):
        return node
    elif isinstance(node, numbers.Integral):
        return LiteralInt(node)
    elif isinstance(node, numbers.Real):
        return LiteralFloat(node)
    else:
        raise RuntimeError("Unexpected LExpr type %s:\n%s" % (type(node), str(node)))


def as_lexpr_or_string_symbol(node):
    if isinstance(node, str):
        return Symbol(node)
    return as_lexpr(node)


def as_symbol(symbol):
    if isinstance(symbol, str):
        symbol = Symbol(symbol)
    assert isinstance(symbol, Symbol)
    return symbol


def flattened_indices(indices, shape):
    """Return a flattened indexing expression.

    Given a tuple of indices and a shape tuple, return
    a CNode expression for flattened indexing into multidimensional
    array.

    Indices and shape entries can be int values, str symbol names, or
    CNode expressions.

    """
    n = len(shape)
    if n == 0:
        # Scalar
        return as_lexpr(0)
    elif n == 1:
        # Simple vector
        return as_lexpr(indices[0])
    else:
        # 2d or higher
        strides = [None] * (n - 2) + [shape[-1], 1]
        for i in range(n - 3, -1, -1):
            strides[i] = Mul(shape[i + 1], strides[i + 1])
        result = indices[-1]
        for i in range(n - 2, -1, -1):
            result = Add(Mul(strides[i], indices[i]), result)
        return result


class Statement(LNode):
    """Make an expression into a statement."""

    is_scoped = False

    def __init__(self, expr):
        self.expr = as_lexpr(expr)

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.expr == other.expr


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

    def __init__(self, typename, symbol, value=None):
        # No type system yet, just using strings
        assert isinstance(typename, str)
        self.typename = typename

        # Allow Symbol or just a string
        self.symbol = as_symbol(symbol)

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

    def __init__(self, typename, symbol, sizes=None, values=None):
        assert isinstance(typename, str)
        self.typename = typename

        if isinstance(symbol, FlattenedArray):
            if sizes is None:
                assert symbol.dims is not None
                sizes = symbol.dims
            elif symbol.dims is not None:
                assert symbol.dims == sizes
            self.symbol = symbol.array
        else:
            self.symbol = as_symbol(symbol)

        if isinstance(sizes, int):
            sizes = (sizes,)
        self.sizes = tuple(sizes)

        # NB! No type checking, assuming nested lists of literal values. Not applying as_lexpr.
        if isinstance(values, (list, tuple)):
            self.values = np.asarray(values)
        else:
            self.values = values

    def __eq__(self, other):
        attributes = ("typename", "symbol", "sizes", "padlen", "values")
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

    def __init__(self, index, begin, end, body, index_type="int"):
        self.index = as_lexpr_or_string_symbol(index)
        self.begin = as_lexpr(begin)
        self.end = as_lexpr(end)
        self.body = as_statement(body)
        self.index_type = index_type

    def __eq__(self, other):
        attributes = ("index", "begin", "end", "body", "index_type")
        return isinstance(other, type(self)) and all(
            getattr(self, name) == getattr(self, name) for name in attributes
        )


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
            "Unexpected CStatement type %s:\n%s" % (type(node), str(node))
        )
