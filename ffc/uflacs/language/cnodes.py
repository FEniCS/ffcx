# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

import numpy
import numbers

from ffc.uflacs.language.format_value import format_value, format_float, format_int
from ffc.uflacs.language.format_lines import format_indented_lines, Indented
from ffc.uflacs.language.precedence import PRECEDENCE


"""CNode TODO:
- Array copy statement
- Extend ArrayDecl and ArrayAccess with support for
  flattened but conceptually multidimensional arrays,
  maybe even with padding (FlattenedArray possibly covers what we need)
- Function declaration
- TypeDef
- Type
- TemplateArgumentList
- Class declaration
- Class definition
"""


############## Some helper functions

def assign_loop(dst, src, ranges):
    """Generate a nested loop over a list of ranges, assigning src to dst in the innermost loop.

    Ranges is a list on the format [(index, begin, end),...].
    """
    code = Assign(dst, src)
    for i, b, e in reversed(ranges):
        code = ForRange(i, b, e, code)
    return code


def accumulate_loop(dst, src, ranges):
    """Generate a nested loop over a list of ranges, adding src to dst in the innermost loop.

    Ranges is a list on the format [(index, begin, end),...].
    """
    code = AssignAdd(dst, src)
    for i, b, e in reversed(ranges):
        code = ForRange(i, b, e, code)
    return code


def scale_loop(dst, factor, ranges):
    """Generate a nested loop over a list of ranges, multiplying dst with factor in the innermost loop.

    Ranges is a list on the format [(index, begin, end),...].
    """
    code = AssignMul(dst, factor)
    for i, b, e in reversed(ranges):
        code = ForRange(i, b, e, code)
    return code


def is_zero_cexpr(cexpr):
    return (
        (isinstance(cexpr, LiteralFloat) and cexpr.value == 0.0)
        or (isinstance(cexpr, LiteralInt) and cexpr.value == 0)
        )


def is_one_cexpr(cexpr):
    return (
        (isinstance(cexpr, LiteralFloat) and cexpr.value == 1.0)
        or (isinstance(cexpr, LiteralInt) and cexpr.value == 1)
        )


def is_negative_one_cexpr(cexpr):
    return (
        (isinstance(cexpr, LiteralFloat) and cexpr.value == -1.0)
        or (isinstance(cexpr, LiteralInt) and cexpr.value == -1)
        )


def float_product(factors):
    "Build product of float factors, simplifying ones and zeros and returning 1.0 if empty sequence."
    factors = [f for f in factors if not is_one_cexpr(f)]
    if len(factors) == 0:
        return LiteralFloat(1.0)
    elif len(factors) == 1:
        return factors[0]
    else:
        for f in factors:
            if is_zero_cexpr(f):
                return f
        return Product(factors)


def MemZeroRange(name, begin, end):
    name = as_cexpr_or_string_symbol(name)
    return Call("std::fill", (name + begin, name + end, LiteralFloat(0.0)))
    #return Call("std::fill", (AddressOf(name[begin]), AddressOf(name[end]), LiteralFloat(0.0)))


def MemZero(name, size):
    name = as_cexpr_or_string_symbol(name)
    size = as_cexpr(size)
    return Call("std::fill_n", (name, size, LiteralFloat(0.0)))


def MemCopy(src, dst, size):
    src = as_cexpr_or_string_symbol(src)
    dst = as_cexpr_or_string_symbol(dst)
    size = as_cexpr(size)
    return Call("std::copy_n", (src, size, dst))


############## CNode core

class CNode(object):
    "Base class for all C AST nodes."
    __slots__ = ()

    def __str__(self):
        name = self.__class__.__name__
        raise NotImplementedError("Missing implementation of __str__ in " + name)

    def __eq__(self, other):
        name = self.__class__.__name__
        raise NotImplementedError("Missing implementation of __eq__ in " + name)

    def __ne__(self, other):
        return not self.__eq__(other)


CNode.debug = False


############## CExpr base classes

class CExpr(CNode):
    """Base class for all C expressions.

    All subtypes should define a 'precedence' class attribute.
    """
    __slots__ = ()

    def ce_format(self, precision=None):
        raise NotImplementedError("Missing implementation of ce_format() in CExpr.")

    def __str__(self):
        try:
            s = self.ce_format()
        except Exception:
            if CNode.debug:
                print("Error in CExpr string formatting. Inspect self.")
                import IPython; IPython.embed()
            raise
        return s

    def __getitem__(self, indices):
        return ArrayAccess(self, indices)

    def __neg__(self):
        if isinstance(self, LiteralFloat):
            return LiteralFloat(-self.value)
        if isinstance(self, LiteralInt):
            return LiteralInt(-self.value)
        return Neg(self)

    def __add__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(self):
            return other
        if is_zero_cexpr(other):
            return self
        if isinstance(other, Neg):
            return Sub(self, other.arg)
        return Add(self, other)

    def __radd__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(self):
            return other
        if is_zero_cexpr(other):
            return self
        if isinstance(self, Neg):
            return Sub(other, self.arg)
        return Add(other, self)

    def __sub__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(self):
            return -other
        if is_zero_cexpr(other):
            return self
        if isinstance(other, Neg):
            return Add(self, other.arg)
        return Sub(self, other)

    def __rsub__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(self):
            return other
        if is_zero_cexpr(other):
            return -self
        if isinstance(self, Neg):
            return Add(other, self.arg)
        return Sub(other, self)

    def __mul__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(self):
            return self
        if is_zero_cexpr(other):
            return other
        if is_one_cexpr(self):
            return other
        if is_one_cexpr(other):
            return self
        if is_negative_one_cexpr(other):
            return Neg(self)
        if is_negative_one_cexpr(self):
            return Neg(other)
        return Mul(self, other)

    def __rmul__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(self):
            return self
        if is_zero_cexpr(other):
            return other
        if is_one_cexpr(self):
            return other
        if is_one_cexpr(other):
            return self
        if is_negative_one_cexpr(other):
            return Neg(self)
        if is_negative_one_cexpr(self):
            return Neg(other)
        return Mul(other, self)

    def __div__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(other):
            raise ValueError("Division by zero!")
        if is_zero_cexpr(self):
            return self
        return Div(self, other)

    def __rdiv__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(self):
            raise ValueError("Division by zero!")
        if is_zero_cexpr(other):
            return other
        return Div(other, self)

    # TODO: Error check types? Can't do that exactly as symbols here have no type.
    __truediv__ = __div__
    __rtruediv__ = __rdiv__
    __floordiv__ = __div__
    __rfloordiv__ = __rdiv__

    def __mod__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(other):
            raise ValueError("Division by zero!")
        if is_zero_cexpr(self):
            return self
        return Mod(self, other)

    def __rmod__(self, other):
        other = as_cexpr(other)
        if is_zero_cexpr(self):
            raise ValueError("Division by zero!")
        if is_zero_cexpr(other):
            return other
        return Mod(other, self)


class CExprOperator(CExpr):
    """Base class for all C expression operator."""
    __slots__ = ()
    sideeffect = False


class CExprTerminal(CExpr):
    """Base class for all C expression terminals."""
    __slots__ = ()
    sideeffect = False


############## CExprTerminal types

class CExprLiteral(CExprTerminal):
    "A float or int literal value."
    __slots__ = ()
    precedence = PRECEDENCE.LITERAL


class Null(CExprLiteral):
    "A null pointer literal."
    __slots__ = ()
    precedence = PRECEDENCE.LITERAL

    def ce_format(self, precision=None):
        # C or old C++ version
        #return "NULL"
        # C++11 version
        return "nullptr"

    def __eq__(self, other):
        return isinstance(other, Null)


class LiteralFloat(CExprLiteral):
    "A floating point literal value."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (float, int, numpy.number))
        self.value = value

    def ce_format(self, precision=None):
        return format_float(self.value, precision)

    def __eq__(self, other):
        return isinstance(other, LiteralFloat) and self.value == other.value

    def __bool__(self):
        return bool(self.value)

    __nonzero__ = __bool__

    def __float__(self):
        return float(self.value)


class LiteralInt(CExprLiteral):
    "An integer literal value."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (int, numpy.number))
        self.value = value

    def ce_format(self, precision=None):
        return str(self.value)

    def __eq__(self, other):
        return isinstance(other, LiteralInt) and self.value == other.value

    def __bool__(self):
        return bool(self.value)

    __nonzero__ = __bool__

    def __int__(self):
        return int(self.value)

    def __float__(self):
        return float(self.value)


class LiteralBool(CExprLiteral):
    "A boolean literal value."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (bool,))
        self.value = value

    def ce_format(self, precision=None):
        return "true" if self.value else "false"

    def __eq__(self, other):
        return isinstance(other, LiteralBool) and self.value == other.value

    def __bool__(self):
        return bool(self.value)

    __nonzero__ = __bool__


class LiteralString(CExprLiteral):
    "A boolean literal value."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (str,))
        assert '"' not in value
        self.value = value

    def ce_format(self, precision=None):
        return '"%s"' % (self.value,)

    def __eq__(self, other):
        return isinstance(other, LiteralString) and self.value == other.value


class Symbol(CExprTerminal):
    "A named symbol."
    __slots__ = ("name",)
    precedence = PRECEDENCE.SYMBOL

    def __init__(self, name):
        assert isinstance(name, str)
        self.name = name

    def ce_format(self, precision=None):
        return self.name

    def __eq__(self, other):
        return isinstance(other, Symbol) and self.name == other.name


class VerbatimExpr(CExprTerminal):
    """A verbatim copy of an expression source string.

    Handled as having the lowest precedence which will introduce parentheses around it most of the time."""
    __slots__ = ("codestring",)
    precedence = PRECEDENCE.LOWEST

    def __init__(self, codestring):
        assert isinstance(codestring, str)
        self.codestring = codestring

    def ce_format(self, precision=None):
        return self.codestring

    def __eq__(self, other):
        return isinstance(other, VerbatimExpr) and self.codestring == other.codestring


class New(CExpr):
    __slots__ = ("typename",)
    def __init__(self, typename):
        assert isinstance(typename, str)
        self.typename = typename

    def ce_format(self, precision=None):
        return "new %s()" % (self.typename,)

    def __eq__(self, other):
        return isinstance(other, New) and self.typename == other.typename


############## CExprOperator base classes

class UnaryOp(CExprOperator):
    "Base class for unary operators."
    __slots__ = ("arg",)
    def __init__(self, arg):
        self.arg = as_cexpr(arg)

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.arg == other.arg


class PrefixUnaryOp(UnaryOp):
    "Base class for prefix unary operators."
    __slots__ = ()
    def ce_format(self, precision=None):
        arg = self.arg.ce_format(precision)
        if self.arg.precedence >= self.precedence:
            arg = '(' + arg + ')'
        return self.op + arg

    def __eq__(self, other):
        return isinstance(other, type(self))


class PostfixUnaryOp(UnaryOp):
    "Base class for postfix unary operators."
    __slots__ = ()
    def ce_format(self, precision=None):
        arg = self.arg.ce_format(precision)
        if self.arg.precedence >= self.precedence:
            arg = '(' + arg + ')'
        return arg + self.op

    def __eq__(self, other):
        return isinstance(other, type(self))


class BinOp(CExprOperator):
    __slots__ = ("lhs", "rhs")
    def __init__(self, lhs, rhs):
        self.lhs = as_cexpr(lhs)
        self.rhs = as_cexpr(rhs)

    def ce_format(self, precision=None):
        # Format children
        lhs = self.lhs.ce_format(precision)
        rhs = self.rhs.ce_format(precision)

        # Apply parentheses
        if self.lhs.precedence > self.precedence:
            lhs = '(' + lhs + ')'
        if self.rhs.precedence >= self.precedence:
            rhs = '(' + rhs + ')'

        # Return combined string
        return lhs + (" " + self.op + " ") + rhs

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.lhs == other.lhs
                    and self.rhs == other.rhs)


class NaryOp(CExprOperator):
    "Base class for special n-ary operators."
    __slots__ = ("args",)
    def __init__(self, args):
        self.args = [as_cexpr(arg) for arg in args]

    def ce_format(self, precision=None):
        # Format children
        args = [arg.ce_format(precision) for arg in self.args]

        # Apply parentheses
        for i in range(len(args)):
            if self.args[i].precedence >= self.precedence:
                args[i] = '(' + args[i] + ')'

        # Return combined string
        op = " " + self.op + " "
        s = args[0]
        for i in range(1, len(args)):
            s += op + args[i]
        return s

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and len(self.args) == len(other.args)
                    and all(a == b for a, b in zip(self.args, other.args)))


############## CExpr unary operators

class Dereference(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.DEREFERENCE
    op = "*"

class AddressOf(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.ADDRESSOF
    op = "&"

class SizeOf(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.SIZEOF
    op = "sizeof"

class Neg(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.NEG
    op = "-"

class Pos(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.POS
    op = "+"

class Not(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.NOT
    op = "!"

class BitNot(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.BIT_NOT
    op = "~"

class PreIncrement(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.PRE_INC
    sideeffect = True
    op = "++"

class PreDecrement(PrefixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.PRE_DEC
    sideeffect = True
    op = "--"

class PostIncrement(PostfixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.POST_INC
    sideeffect = True
    op = "++"

class PostDecrement(PostfixUnaryOp):
    __slots__ = ()
    precedence = PRECEDENCE.POST_DEC
    sideeffect = True
    op = "--"

############## CExpr binary operators

class Add(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.ADD
    op = "+"

class Sub(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.SUB
    op = "-"

class Mul(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.MUL
    op = "*"

class Div(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.DIV
    op = "/"

class Mod(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.MOD
    op = "%"

class EQ(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.EQ
    op = "=="

class NE(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.NE
    op = "!="

class LT(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.LT
    op = "<"

class GT(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.GT
    op = ">"

class LE(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.LE
    op = "<="

class GE(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.GE
    op = ">="

class And(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.AND
    op = "&&"

class Or(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.OR
    op = "||"

class BitAnd(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.BIT_AND
    op = "&"

class BitXor(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.BIT_XOR
    op = "^"

class BitOr(BinOp):
    __slots__ = ()
    precedence = PRECEDENCE.BIT_OR
    op = "|"

class Sum(NaryOp):
    "Sum of any number of operands."
    __slots__ = ()
    precedence = PRECEDENCE.ADD
    op = "+"

class Product(NaryOp):
    "Product of any number of operands."
    __slots__ = ()
    precedence = PRECEDENCE.MUL
    op = "*"

class AssignOp(BinOp):
    "Base class for assignment operators."
    __slots__ = ()
    precedence = PRECEDENCE.ASSIGN
    sideeffect = True
    def __init__(self, lhs, rhs):
        BinOp.__init__(self, as_cexpr_or_string_symbol(lhs), rhs)

class Assign(AssignOp):
    __slots__ = ()
    op = "="

class AssignAdd(AssignOp):
    __slots__ = ()
    op = "+="

class AssignSub(AssignOp):
    __slots__ = ()
    op = "-="

class AssignMul(AssignOp):
    __slots__ = ()
    op = "*="

class AssignDiv(AssignOp):
    __slots__ = ()
    op = "/="

class AssignMod(AssignOp):
    __slots__ = ()
    op = "%="

class AssignLShift(AssignOp):
    __slots__ = ()
    op = "<<="

class AssignRShift(AssignOp):
    __slots__ = ()
    op = ">>="

class AssignAnd(AssignOp):
    __slots__ = ()
    op = "&&="

class AssignOr(AssignOp):
    __slots__ = ()
    op = "||="

class AssignBitAnd(AssignOp):
    __slots__ = ()
    op = "&="

class AssignBitXor(AssignOp):
    __slots__ = ()
    op = "^="

class AssignBitOr(AssignOp):
    __slots__ = ()
    op = "|="


############## CExpr operators

class FlattenedArray(object):
    """Syntax carrying object only, will get translated on __getitem__ to ArrayAccess."""
    __slots__ = ("array", "strides", "offset", "dims")
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
            dims = tuple(as_cexpr(i) for i in dims)
            self.dims = dims
            n = len(dims)
            literal_one = LiteralInt(1)
            strides = [literal_one]*n
            for i in range(n-2, -1, -1):
                s = strides[i+1]
                d = dims[i+1]
                if d == literal_one:
                    strides[i] = s
                elif s == literal_one:
                    strides[i] = d
                else:
                    strides[i] = d * s
        else:
            self.dims = None
            assert isinstance(strides, (list, tuple))
            strides = tuple(as_cexpr(i) for i in strides)
        self.strides = strides
        self.offset = None if offset is None else as_cexpr(offset)

    def __getitem__(self, indices):
        if not isinstance(indices, (list,tuple)):
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
            flat = (i if s == literal_one else s * i)
            if self.offset is not None:
                flat = self.offset + flat
            for i, s in zip(indices[1:n], self.strides[1:n]):
                flat = flat + (i if s == literal_one else s * i)
        # Delay applying ArrayAccess until we have all indices
        if n == len(self.strides):
            return ArrayAccess(self.array, flat)
        else:
            return FlattenedArray(self.array, strides=self.strides[n:], offset=flat)


class ArrayAccess(CExprOperator):
    __slots__ = ("array", "indices")
    precedence = PRECEDENCE.SUBSCRIPT

    def __init__(self, array, indices):
        # Typecheck array argument
        if isinstance(array, str):
            array = Symbol(array)
        if isinstance(array, Symbol):
            self.array = array
        elif isinstance(array, ArrayDecl):
            self.array = array.symbol
        else:
            raise ValueError("Unexpected array type %s." % (type(array).__name__,))

        # Allow expressions or literals as indices
        if not isinstance(indices, (list, tuple)):
            indices = (indices,)
        self.indices = tuple(as_cexpr_or_string_symbol(i) for i in indices)

        # Early error checking for negative array dimensions
        if any(isinstance(i, int) and i < 0 for i in self.indices):
            raise ValueError("Index value < 0.")

        # Additional dimension checks possible if we get an ArrayDecl instead of just a name
        if isinstance(array, ArrayDecl):
            if len(self.indices) != len(array.sizes):
                raise ValueError("Invalid number of indices.")
            ints = (int, LiteralInt)
            if any((isinstance(i, ints) and isinstance(d, ints) and int(i) >= int(d))
                   for i, d in zip(self.indices, array.sizes)):
                raise ValueError("Index value >= array dimension.")

    def __getitem__(self, indices):
        "Handling nested expr[i][j]."
        if isinstance(indices, list):
            indices = tuple(indices)
        elif not isinstance(indices, tuple):
            indices = (indices,)
        return ArrayAccess(self.array, self.indices + indices)

    def ce_format(self, precision=None):
        s = self.array.ce_format(precision)
        for index in self.indices:
            s += "[" + index.ce_format(precision) + "]"
        return s

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.array == other.array
                    and self.indices == other.indices)


class Conditional(CExprOperator):
    __slots__ = ("condition", "true", "false")
    precedence = PRECEDENCE.CONDITIONAL

    def __init__(self, condition, true, false):
        self.condition = as_cexpr(condition)
        self.true = as_cexpr(true)
        self.false = as_cexpr(false)

    def ce_format(self, precision=None):
        # Format children
        c = self.condition.ce_format(precision)
        t = self.true.ce_format(precision)
        f = self.false.ce_format(precision)

        # Apply parentheses
        if self.condition.precedence >= self.precedence:
            c = '(' + c + ')'
        if self.true.precedence >= self.precedence:
            t = '(' + t + ')'
        if self.false.precedence >= self.precedence:
            f = '(' + f + ')'

        # Return combined string
        return c + " ? " + t + " : " + f

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.condition == other.condition
                    and self.true == other.true
                    and self.false == other.false)


class Call(CExprOperator):
    __slots__ = ("function", "arguments")
    precedence = PRECEDENCE.CALL
    sideeffect = True

    def __init__(self, function, arguments=None):
        self.function = as_cexpr_or_string_symbol(function)

        # Accept None, single, or multple arguments; literals or CExprs
        if arguments is None:
            arguments = ()
        elif not isinstance(arguments, (tuple, list)):
            arguments = (arguments,)
        self.arguments = [as_cexpr(arg) for arg in arguments]

    def ce_format(self, precision=None):
        args = ", ".join(arg.ce_format(precision) for arg in self.arguments)
        return self.function.ce_format(precision) + "(" + args + ")"

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.function == other.function
                    and self.arguments == other.arguments)

def Sqrt(x):
    return Call("std::sqrt", x)


############## Convertion function to expression nodes

def _is_zero_valued(values):
    if isinstance(values, (numbers.Integral, LiteralInt)):
        return int(values) == 0
    elif isinstance(values, (numbers.Number, LiteralFloat)):
        return float(values) == 0.0
    else:
        return numpy.count_nonzero(values) == 0


def as_cexpr(node):
    """Typechecks and wraps an object as a valid CExpr.

    Accepts CExpr nodes, treats int and float as literals, and treats a string as a symbol.
    """
    if isinstance(node, CExpr):
        return node
    elif isinstance(node, bool):
        return LiteralBool(node)
    elif isinstance(node, numbers.Integral):
        return LiteralInt(node)
    elif isinstance(node, numbers.Real):
        return LiteralFloat(node)
    elif isinstance(node, str):
        raise RuntimeError("Got string for CExpr, this is ambiguous: %s" % (node,))
    else:
        raise RuntimeError("Unexpected CExpr type %s:\n%s" % (type(node), str(node)))


def as_cexpr_or_string_symbol(node):
    if isinstance(node, str):
        return Symbol(node)
    return as_cexpr(node)


def as_cexpr_or_verbatim(node):
    if isinstance(node, str):
        return VerbatimExpr(node)
    return as_cexpr(node)


def as_cexpr_or_literal(node):
    if isinstance(node, str):
        return LiteralString(node)
    return as_cexpr(node)


def as_symbol(symbol):
    if isinstance(symbol, str):
        symbol = Symbol(symbol)
    assert isinstance(symbol, Symbol)
    return symbol


def flattened_indices(indices, shape):
    """Given a tuple of indices and a shape tuple,
    return CNode expression for flattened indexing
    into multidimensional array.

    Indices and shape entries can be int values, str symbol names, or CNode expressions.
    """
    n = len(shape)
    if n == 0:
        # Scalar
        return as_cexpr(0)
    elif n == 1:
        # Simple vector
        return as_cexpr(indices[0])
    else:
        # 2d or higher
        strides = [None]*(n-2) + [shape[-1], 1]
        for i in range(n-3, -1, -1):
            strides[i] = Mul(shape[i+1], strides[i+1])
        result = indices[-1]
        for i in range(n-2, -1, -1):
            result = Add(Mul(strides[i], indices[i]), result)
        return result


############## Base class for all statements

class CStatement(CNode):
    """Base class for all C statements.

    Subtypes do _not_ define a 'precedence' class attribute.
    """
    __slots__ = ()

    # True if statement contains its own scope, false by default to be on the safe side
    is_scoped = False

    def cs_format(self, precision=None):
        "Return S: string | list(S) | Indented(S)."
        raise NotImplementedError("Missing implementation of cs_format() in CStatement.")

    def __str__(self):
        try:
            s = self.cs_format()
        except Exception:
            if CNode.debug:
                print("Error in CStatement string formatting. Inspect self.")
                import IPython; IPython.embed()
            raise
        return format_indented_lines(s)


############## Statements

class VerbatimStatement(CStatement):
    "Wraps a source code string to be pasted verbatim into the source code."
    __slots__ = ("codestring",)
    is_scoped = False

    def __init__(self, codestring):
        assert isinstance(codestring, str)
        self.codestring = codestring

    def cs_format(self, precision=None):
        return self.codestring

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.codestring == other.codestring)


class Statement(CStatement):
    "Make an expression into a statement."
    __slots__ = ("expr",)
    is_scoped = False
    def __init__(self, expr):
        self.expr = as_cexpr(expr)

    def cs_format(self, precision=None):
        return self.expr.ce_format(precision) + ";"

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.expr == other.expr)


class StatementList(CStatement):
    "A simple sequence of statements. No new scopes are introduced."
    __slots__ = ("statements",)

    def __init__(self, statements):
        self.statements = [as_cstatement(st) for st in statements]

    @property
    def is_scoped(self):
        return all(st.is_scoped for st in self.statements)

    def cs_format(self, precision=None):
        return [st.cs_format(precision) for st in self.statements]

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.statements == other.statements)


############## Simple statements

class Using(CStatement):
    __slots__ = ("name",)
    is_scoped = True
    def __init__(self, name):
        assert isinstance(name, str)
        self.name = name

    def cs_format(self, precision=None):
        return "using " + self.name + ";"

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.name == other.name)


class Break(CStatement):
    __slots__ = ()
    is_scoped = True
    def cs_format(self, precision=None):
        return "break;"

    def __eq__(self, other):
        return isinstance(other, type(self))


class Continue(CStatement):
    __slots__ = ()
    is_scoped = True
    def cs_format(self, precision=None):
        return "continue;"

    def __eq__(self, other):
        return isinstance(other, type(self))


class Return(CStatement):
    __slots__ = ("value",)
    is_scoped = True
    def __init__(self, value=None):
        if value is None:
            self.value = None
        else:
            self.value = as_cexpr(value)

    def cs_format(self, precision=None):
        if self.value is None:
            return "return;"
        else:
            return "return %s;" % (self.value.ce_format(precision),)

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.value == other.value)


class Case(CStatement):
    __slots__ = ("value",)
    is_scoped = False
    def __init__(self, value):
        # NB! This is too permissive and will allow invalid case arguments.
        self.value = as_cexpr(value)

    def cs_format(self, precision=None):
        return "case " + self.value.ce_format(precision) + ":"

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.value == other.value)


class Default(CStatement):
    __slots__ = ()
    is_scoped = False
    def cs_format(self, precision=None):
        return "default:"

    def __eq__(self, other):
        return isinstance(other, type(self))


class Throw(CStatement):
    __slots__ = ("exception", "message")
    is_scoped = True
    def __init__(self, exception, message):
        assert isinstance(exception, str)
        assert isinstance(message, str)
        self.exception = exception
        self.message = message

    def cs_format(self, precision=None):
        assert '"' not in self.message
        return "throw " + self.exception + '("' + self.message + '");'

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.message == other.message
                    and self.exception == other.exception)


class Comment(CStatement):
    "Line comment(s) used for annotating the generated code with human readable remarks."
    __slots__ = ("comment",)
    is_scoped = True
    def __init__(self, comment):
        assert isinstance(comment, str)
        self.comment = comment

    def cs_format(self, precision=None):
        lines = self.comment.strip().split("\n")
        return ["// " + line.strip() for line in lines]

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.comment == other.comment)


def NoOp():
    return Comment("Do nothing")


def commented_code_list(code, comments):
    "Convenience wrapper for adding comment to code list if the list is not empty."
    if isinstance(code, CNode):
        code = [code]
    assert isinstance(code, list)
    if code:
        if not isinstance(comments, (list, tuple)):
            comments = [comments]
        comments = [Comment(c) for c in comments]
        code = comments + code
    return code


class Pragma(CStatement):
    "Pragma comments used for compiler-specific annotations."
    __slots__ = ("comment",)
    is_scoped = True
    def __init__(self, comment):
        assert isinstance(comment, str)
        self.comment = comment

    def cs_format(self, precision=None):
        assert "\n" not in self.comment
        return "#pragma " + self.comment

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.comment == other.comment)


############## Type and variable declarations

class VariableDecl(CStatement):
    "Declare a variable, optionally define initial value."
    __slots__ = ("typename", "symbol", "value")
    is_scoped = False
    def __init__(self, typename, symbol, value=None):

        # No type system yet, just using strings
        assert isinstance(typename, str)
        self.typename = typename

        # Allow Symbol or just a string
        self.symbol = as_symbol(symbol)

        if value is not None:
            value = as_cexpr(value)
        self.value = value

    def cs_format(self, precision=None):
        code = self.typename + " " + self.symbol.name
        if self.value is not None:
            code += " = " + self.value.ce_format(precision)
        return code + ";"

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.typename == other.typename
                    and self.symbol == other.symbol
                    and self.value == other.value)


def leftover(size, padlen):
    "Return minimum integer to add to size to make it divisible by padlen."
    return (padlen - (size % padlen)) % padlen


def pad_dim(dim, padlen):
    "Make dim divisible by padlen."
    return ((dim + padlen - 1) // padlen) * padlen


def pad_innermost_dim(shape, padlen):
    "Make the last dimension in shape divisible by padlen."
    if not shape:
        return ()
    shape = list(shape)
    if padlen:
        shape[-1] = pad_dim(shape[-1], padlen)
    return tuple(shape)


def build_1d_initializer_list(values, formatter, padlen=0, precision=None):
    '''Return a list containing a single line formatted like "{ 0.0, 1.0, 2.0 }"'''
    if formatter == str:
        formatter = lambda x, p: str(x)
    tokens = ["{ "]
    if numpy.product(values.shape) > 0:
        sep = ", "
        fvalues = [formatter(v, precision) for v in values]
        for v in fvalues[:-1]:
            tokens.append(v)
            tokens.append(sep)
        tokens.append(fvalues[-1])
        if padlen:
            # Add padding
            zero = formatter(values.dtype.type(0), precision)
            for i in range(leftover(len(values), padlen)):
                tokens.append(sep)
                tokens.append(zero)
    tokens += " }"
    return "".join(tokens)


def build_initializer_lists(values, sizes, level, formatter, padlen=0, precision=None):
    """Return a list of lines with initializer lists for a multidimensional array.

    Example output::

        { { 0.0, 0.1 },
          { 1.0, 1.1 } }
    """
    if formatter == str:
        formatter = lambda x, p: str(x)
    values = numpy.asarray(values)
    assert numpy.product(values.shape) == numpy.product(sizes)
    assert len(sizes) > 0
    assert len(values.shape) > 0
    assert len(sizes) == len(values.shape)
    assert numpy.all(values.shape == sizes)

    r = len(sizes)
    assert r > 0
    if r == 1:
        return [build_1d_initializer_list(values, formatter, padlen=padlen, precision=precision)]
    else:
        # Render all sublists
        parts = []
        for val in values:
            sublist = build_initializer_lists(val, sizes[1:], level+1, formatter, padlen=padlen, precision=precision)
            parts.append(sublist)
        # Add comma after last line in each part except the last one
        for part in parts[:-1]:
            part[-1] += ","
        # Collect all lines in flat list
        lines = []
        for part in parts:
            lines.extend(part)
        # Enclose lines in '{ ' and ' }' and indent lines in between
        lines[0] = "{ " + lines[0]
        for i in range(1,len(lines)):
            lines[i] = "  " + lines[i]
        lines[-1] += " }"
        return lines


class ArrayDecl(CStatement):
    """A declaration or definition of an array.

    Note that just setting values=0 is sufficient
    to initialize the entire array to zero.

    Otherwise use nested lists of lists to represent
    multidimensional array values to initialize to.
    """
    __slots__ = ("typename", "symbol", "sizes", "alignas", "padlen", "values")
    is_scoped = False
    def __init__(self, typename, symbol, sizes=None, values=None, alignas=None, padlen=0):
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

        # NB! No type checking, assuming nested lists of literal values. Not applying as_cexpr.
        if isinstance(values, (list, tuple)):
            self.values = numpy.asarray(values)
        else:
            self.values = values

        self.alignas = alignas
        self.padlen = padlen

    def __getitem__(self, indices):
        """Allow using array declaration object as the array when indexed.

        A = ArrayDecl("int", "A", (2,3))
        code = [A, Assign(A[0,0], 1.0)]
        """
        return ArrayAccess(self, indices)

    def cs_format(self, precision=None):
        # Pad innermost array dimension
        sizes = pad_innermost_dim(self.sizes, self.padlen)

        # Add brackets
        brackets = ''.join("[%d]" % n for n in sizes)

        # Join declaration
        decl = self.typename + " " + self.symbol.name + brackets

        # NB! C++11 style alignas prefix syntax.
        # If trying other target languages, must use other syntax.
        if self.alignas:
            align = "alignas(%d)" % int(self.alignas)
            decl = align + " " + decl

        if self.values is None:
            # Undefined initial values
            return decl + ";"
        elif _is_zero_valued(self.values):
            # Zero initial values
            # (NB! C++ style zero initialization, not sure about other target languages)
            return decl + " = {};"
        else:
            # Construct initializer lists for arbitrary multidimensional array values
            if self.values.dtype.kind == "f":
                formatter = format_float
            elif self.values.dtype.kind == "i":
                formatter = format_int
            else:
                formatter = format_value
            initializer_lists = build_initializer_lists(self.values, self.sizes, 0,
                                                        formatter, padlen=self.padlen,
                                                        precision=precision)
            if len(initializer_lists) == 1:
                return decl + " = " + initializer_lists[0] + ";"
            else:
                initializer_lists[-1] += ";" # Close statement on final line
                return (decl + " =", Indented(initializer_lists))

    def __eq__(self, other):
        attributes = ("typename", "symbol", "sizes", "alignas", "padlen", "values")
        return (isinstance(other, type(self))
                    and all(getattr(self, name) == getattr(self, name)
                            for name in attributes))



############## Scoped statements

class Scope(CStatement):
    __slots__ = ("body",)
    is_scoped = True
    def __init__(self, body):
        self.body = as_cstatement(body)

    def cs_format(self, precision=None):
        return ("{", Indented(self.body.cs_format(precision)), "}")

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.body == other.body)


class Namespace(CStatement):
    __slots__ = ("name", "body")
    is_scoped = True
    def __init__(self, name, body):
        assert isinstance(name, str)
        self.name = name
        self.body = as_cstatement(body)

    def cs_format(self, precision=None):
        return ("namespace " + self.name,
                "{", Indented(self.body.cs_format(precision)), "}")

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.name == other.name
                    and self.body == other.body)


def _is_scoped_statement(body):
    return


def _is_simple_if_body(body):
    if isinstance(body, StatementList):
        if len(body.statements) > 1:
            return False
        body, = body.statements
    return isinstance(body, (Return, AssignOp, Break, Continue))


class If(CStatement):
    __slots__ = ("condition", "body")
    is_scoped = True
    def __init__(self, condition, body):
        self.condition = as_cexpr(condition)
        self.body = as_cstatement(body)

    def cs_format(self, precision=None):
        statement = "if (" + self.condition.ce_format(precision) + ")"
        body_fmt = Indented(self.body.cs_format(precision))
        if _is_simple_if_body(self.body):
            return (statement, body_fmt)
        else:
            return (statement, "{", body_fmt, "}")

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.condition == other.condition
                    and self.body == other.body)


class ElseIf(CStatement):
    __slots__ = ("condition", "body")
    is_scoped = True
    def __init__(self, condition, body):
        self.condition = as_cexpr(condition)
        self.body = as_cstatement(body)

    def cs_format(self, precision=None):
        statement = "else if (" + self.condition.ce_format(precision) + ")"
        body_fmt = Indented(self.body.cs_format(precision))
        if _is_simple_if_body(self.body):
            return (statement, body_fmt)
        else:
            return (statement, "{", body_fmt, "}")

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.condition == other.condition
                    and self.body == other.body)


class Else(CStatement):
    __slots__ = ("body",)
    is_scoped = True
    def __init__(self, body):
        self.body = as_cstatement(body)

    def cs_format(self, precision=None):
        statement = "else"
        body_fmt = Indented(self.body.cs_format(precision))
        if _is_simple_if_body(self.body):
            return (statement, body_fmt)
        else:
            return (statement, "{", body_fmt, "}")

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.body == other.body)


class While(CStatement):
    __slots__ = ("condition", "body")
    is_scoped = True
    def __init__(self, condition, body):
        self.condition = as_cexpr(condition)
        self.body = as_cstatement(body)

    def cs_format(self, precision=None):
        return ("while (" + self.condition.ce_format(precision) + ")",
                "{", Indented(self.body.cs_format(precision)), "}")

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.condition == other.condition
                    and self.body == other.body)


class Do(CStatement):
    __slots__ = ("condition", "body")
    is_scoped = True
    def __init__(self, condition, body):
        self.condition = as_cexpr(condition)
        self.body = as_cstatement(body)

    def cs_format(self, precision=None):
        return ("do", "{", Indented(self.body.cs_format(precision)),
                "} while (" + self.condition.ce_format(precision) + ");")

    def __eq__(self, other):
        return (isinstance(other, type(self))
                    and self.condition == other.condition
                    and self.body == other.body)


def as_pragma(pragma):
    if isinstance(pragma, str):
        return Pragma(pragma)
    elif isinstance(pragma, Pragma):
        return pragma
    return None


def is_simple_inner_loop(code):
    if isinstance(code, (ForRange, For)) and code.pragma is None and is_simple_inner_loop(code.body):
        return True
    if isinstance(code, Statement) and isinstance(code.expr, AssignOp):
        return True
    return False


class For(CStatement):
    __slots__ = ("init", "check", "update", "body", "pragma")
    is_scoped = True
    def __init__(self, init, check, update, body, pragma=None):
        self.init = as_cstatement(init)
        self.check = as_cexpr_or_verbatim(check)
        self.update = as_cexpr_or_verbatim(update)
        self.body = as_cstatement(body)
        self.pragma = as_pragma(pragma)

    def cs_format(self, precision=None):
        # The C model here is a bit crude and this causes trouble
        # in the init statement/expression here:
        init = self.init.cs_format(precision)
        assert isinstance(init, str)
        init = init.rstrip(" ;")

        check = self.check.ce_format(precision)
        update = self.update.ce_format(precision)

        prelude = "for (" + init + "; " + check + "; " + update + ")"
        body = Indented(self.body.cs_format(precision))

        # Reduce size of code with lots of simple loops by dropping {} in obviously safe cases
        if is_simple_inner_loop(self.body):
            code = (prelude, body)
        else:
            code = (prelude, "{", body, "}")

        # Add pragma prefix if requested
        if self.pragma is not None:
            code = (self.pragma.cs_format(),) + code

        return code

    def __eq__(self, other):
        attributes = ("init", "check", "update", "body")
        return (isinstance(other, type(self))
                    and all(getattr(self, name) == getattr(self, name)
                            for name in attributes))


class Switch(CStatement):
    __slots__ = ("arg", "cases", "default", "autobreak", "autoscope")
    is_scoped = True
    def __init__(self, arg, cases, default=None, autobreak=True, autoscope=True):
        self.arg = as_cexpr_or_string_symbol(arg)
        self.cases = [(as_cexpr(value), as_cstatement(body)) for value, body in cases]
        if default is not None:
            default = as_cstatement(default)
            defcase = [(None, default)]
        else:
            defcase = []
        self.default = default
        # If this is a switch where every case returns, scopes or breaks are never needed
        if all(isinstance(case[1], Return) for case in self.cases + defcase):
            autobreak = False
            autoscope = False
        if all(case[1].is_scoped for case in self.cases + defcase):
            autoscope = False
        assert autobreak in (True, False)
        assert autoscope in (True, False)
        self.autobreak = autobreak
        self.autoscope = autoscope

    def cs_format(self, precision=None):
        cases = []
        for case in self.cases:
            caseheader = "case " + case[0].ce_format(precision) + ":"
            casebody = case[1].cs_format(precision)
            if self.autoscope:
                casebody = ("{", Indented(casebody), "}")
            if self.autobreak:
                casebody = (casebody, "break;")
            cases.extend([caseheader, Indented(casebody)])

        if self.default is not None:
            caseheader = "default:"
            casebody = self.default.cs_format(precision)
            if self.autoscope:
                casebody = ("{", Indented(casebody), "}")
            cases.extend([caseheader, Indented(casebody)])

        return ("switch (" + self.arg.ce_format(precision) + ")",
                "{", cases, "}")

    def __eq__(self, other):
        attributes = ("arg", "cases", "default", "autobreak", "autoscope")
        return (isinstance(other, type(self))
                    and all(getattr(self, name) == getattr(self, name)
                            for name in attributes))


class ForRange(CStatement):
    "Slightly higher-level for loop assuming incrementing an index over a range."
    __slots__ = ("index", "begin", "end", "body", "pragma", "index_type")
    is_scoped = True
    def __init__(self, index, begin, end, body, index_type="int", vectorize=None):
        self.index = as_cexpr_or_string_symbol(index)
        self.begin = as_cexpr(begin)
        self.end = as_cexpr(end)
        self.body = as_cstatement(body)

        if vectorize:
            pragma = Pragma("omp simd")
        else:
            pragma = None
        self.pragma = pragma

        self.index_type = index_type

    def cs_format(self, precision=None):
        indextype = self.index_type
        index = self.index.ce_format(precision)
        begin = self.begin.ce_format(precision)
        end = self.end.ce_format(precision)

        init = indextype + " " + index + " = " + begin
        check = index + " < " + end
        update = "++" + index

        prelude = "for (" + init + "; " + check + "; " + update + ")"
        body = Indented(self.body.cs_format(precision))

        # Reduce size of code with lots of simple loops by dropping {} in obviously safe cases
        if is_simple_inner_loop(self.body):
            code = (prelude, body)
        else:
            code = (prelude, "{", body, "}")

        # Add vectorization hint if requested
        if self.pragma is not None:
            code = (self.pragma.cs_format(),) + code

        return code

    def __eq__(self, other):
        attributes = ("index", "begin", "end", "body", "pragma", "index_type")
        return (isinstance(other, type(self))
                    and all(getattr(self, name) == getattr(self, name)
                            for name in attributes))


def ForRanges(*ranges, **kwargs):
    ranges = list(reversed(ranges))
    code = kwargs["body"]
    for r in ranges:
        kwargs["body"] = code
        code = ForRange(*r, **kwargs)
    return code


############## Convertion function to statement nodes

def as_cstatement(node):
    "Perform type checking on node and wrap in a suitable statement type if necessary."
    if isinstance(node, StatementList) and len(node.statements) == 1:
        # Cleans up the expression tree a bit
        return node.statements[0]
    elif isinstance(node, CStatement):
        # No-op
        return node
    elif isinstance(node, CExprOperator):
        if node.sideeffect:
            # Special case for using assignment expressions as statements
            return Statement(node)
        else:
            raise RuntimeError(
                "Trying to create a statement of CExprOperator type %s:\n%s"
                % (type(node), str(node)))
    elif isinstance(node, list):
        # Convenience case for list of statements
        if len(node) == 1:
            # Cleans up the expression tree a bit
            return as_cstatement(node[0])
        else:
            return StatementList(node)
    elif isinstance(node, str):
        # Backdoor for flexibility in code generation to allow verbatim pasted statements
        return VerbatimStatement(node)
    else:
        raise RuntimeError("Unexpected CStatement type %s:\n%s" % (type(node), str(node)))
