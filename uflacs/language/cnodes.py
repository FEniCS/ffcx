# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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

from __future__ import print_function # used in some debugging

import numpy

from uflacs.language.format_value import format_value, format_float
from uflacs.language.format_lines import format_indented_lines, Indented
from uflacs.language.precedence import PRECEDENCE

"""CNode TODO:
- Array copy statement
- Memzero statement
- Extend ArrayDecl and ArrayAccess with support for
  flattened but conceptually multidimensional arrays,
  maybe even with padding (FlattenedArray possibly covers what we need)
- ArrayDecl using std::array
- Function declaration
- TypeDef
- Type
- TemplateArgumentList
- Class declaration
- Class definition
"""


############## Some helper functions

def assign_loop(src, dst, ranges):
    """Generate a nested loop over a list of ranges, assigning dst to src in the innermost loop.

    Ranges is a list on the format [(index, begin, end),...].
    """
    code = Assign(src, dst)
    for i, b, e in reversed(ranges):
        code = ForRange(i, b, e, code)
    return code

def accumulate_loop(src, dst, ranges):
    """Generate a nested loop over a list of ranges, adding dst to src in the innermost loop.

    Ranges is a list on the format [(index, begin, end),...].
    """
    code = AssignAdd(src, dst)
    for i, b, e in reversed(ranges):
        code = ForRange(i, b, e, code)
    return code

def scale_loop(src, dst, ranges):
    """Generate a nested loop over a list of ranges, multiplying dst with src in the innermost loop.

    Ranges is a list on the format [(index, begin, end),...].
    """
    code = AssignMul(src, dst)
    for i, b, e in reversed(ranges):
        code = ForRange(i, b, e, code)
    return code


############## CNode core

class CNode(object):
    "Base class for all C AST nodes."
    __slots__ = ()

    def __str__(self):
        name = self.__class__.__name__
        raise NotImplementedError("Missing implementation of __str__ in " + name)

CNode.debug = False


############## CExpr base classes

class CExpr(CNode):
    """Base class for all C expressions.

    All subtypes should define a 'precedence' class attribute.
    """
    __slots__ = ()

    def ce_format(self):
        raise NotImplementedError("Missing implementation of ce_format() in CExpr.")

    def __str__(self):
        try:
            s = self.ce_format()
        except:
            if CNode.debug:
                print("Error in CExpr string formatting. Inspect self.")
                import IPython; IPython.embed()
            raise
        return s

    def __getitem__(self, indices):
        return ArrayAccess(self, indices)

    def __neg__(self):
        return Neg(self)

    def __add__(self, other):
        return Add(self, other)

    def __radd__(self, other):
        return Add(other, self)

    def __sub__(self, other):
        return Sub(self, other)

    def __rsub__(self, other):
        return Sub(other, self)

    def __mul__(self, other):
        return Mul(self, other)

    def __rmul__(self, other):
        return Mul(other, self)

    def __div__(self, other):
        return Div(self, other)

    def __rdiv__(self, other):
        return Div(other, self)

    def __truediv__(self, other):
        return Div(self, other)

    def __rtruediv__(self, other):
        return Div(other, self)

    def __floordiv__(self, other):
        return NotImplemented

    def __rfloordiv__(self, other):
        return NotImplemented

class CExprOperator(CExpr):
    """Base class for all C expression operator."""
    __slots__ = ("children",)
    sideeffect = False

class CExprTerminal(CExpr):
    """Base class for all C expression terminals."""
    __slots__ = ()
    children = ()


############## CExprTerminal types

class CExprLiteral(CExprTerminal):
    "A float or int literal value."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

class Null(CExprLiteral):
    "A null pointer literal."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

    def ce_format(self):
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

    def ce_format(self):
        return format_float(self.value)

    def __eq__(self, other):
        return isinstance(other, LiteralFloat) and self.value == other.value

    def __nonzero__(self):
        return bool(self.value)

    def __float__(self):
        return float(self.value)

class LiteralInt(CExprLiteral):
    "An integer literal value."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (int, numpy.number))
        self.value = value

    def ce_format(self):
        return str(self.value)

    def __eq__(self, other):
        return isinstance(other, LiteralInt) and self.value == other.value

    def __nonzero__(self):
        return bool(self.value)

    def __int__(self):
        return int(self.value)

class LiteralBool(CExprLiteral):
    "A boolean literal value."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (bool,))
        self.value = value

    def ce_format(self):
        return "true" if self.value else "false"

    def __eq__(self, other):
        return isinstance(other, LiteralBool) and self.value == other.value

    def __nonzero__(self):
        return bool(self.value)

    def __bool__(self):
        return bool(self.value)

class LiteralString(CExprLiteral):
    "A boolean literal value."
    __slots__ = ("value",)
    precedence = PRECEDENCE.LITERAL

    def __init__(self, value):
        assert isinstance(value, (str,))
        assert '"' not in value
        self.value = value

    def ce_format(self):
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

    def ce_format(self):
        return self.name

class VerbatimExpr(CExprTerminal):
    """A verbatim copy of an expression source string.

    Handled as having the lowest precedence which will introduce parentheses around it most of the time."""
    __slots__ = ("codestring",)
    precedence = PRECEDENCE.LOWEST

    def __init__(self, codestring):
        assert isinstance(codestring, str)
        self.codestring = codestring

    def ce_format(self):
        return self.codestring

class New(CExpr):
    __slots__ = ("typename",)
    def __init__(self, typename):
        assert isinstance(typename, str)
        self.typename = typename

    def ce_format(self):
        return "new %s()" % (self.typename,)


############## CExprOperator base classes

class UnaryOp(CExprOperator):
    "Base class for unary operators."
    __slots__ = ("arg",)
    def __init__(self, arg):
        self.arg = as_cexpr(arg)

class PrefixUnaryOp(UnaryOp):
    "Base class for prefix unary operators."
    __slots__ = ()
    def ce_format(self):
        arg = self.arg.ce_format()
        if self.arg.precedence >= self.precedence:
            arg = '(' + arg + ')'
        return self.op + arg

class PostfixUnaryOp(UnaryOp):
    "Base class for postfix unary operators."
    __slots__ = ()
    def ce_format(self):
        arg = self.arg.ce_format()
        if self.arg.precedence >= self.precedence:
            arg = '(' + arg + ')'
        return arg + self.op

class BinOp(CExprOperator):
    __slots__ = ("lhs", "rhs")
    def __init__(self, lhs, rhs):
        self.lhs = as_cexpr(lhs)
        self.rhs = as_cexpr(rhs)

    def ce_format(self):
        # Format children
        lhs = self.lhs.ce_format()
        rhs = self.rhs.ce_format()

        # Apply parentheses
        if self.lhs.precedence > self.precedence:
            lhs = '(' + lhs + ')'
        if self.rhs.precedence >= self.precedence:
            rhs = '(' + rhs + ')'

        # Return combined string
        return lhs + (" " + self.op + " ") + rhs

class NaryOp(CExprOperator):
    "Base class for special n-ary operators."
    __slots__ = ("args",)
    def __init__(self, args):
        self.args = [as_cexpr(arg) for arg in args]

    def ce_format(self):
        # Format children
        args = [arg.ce_format() for arg in self.args]

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
    __slots__ = ("array", "strides", "offset")
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
            assert isinstance(strides, (list, tuple))
            strides = tuple(as_cexpr(i) for i in strides)
        self.strides = strides
        self.offset = None if offset is None else as_cexpr(offset)

    def __getitem__(self, indices):
        if not isinstance(indices, (list,tuple)):
            indices = (indices,)
        n = len(indices)
        i, s = (indices[0], self.strides[0])
        literal_one = LiteralInt(1)

        flat = (i if s == literal_one else s * i)
        if self.offset is not None:
            flat = self.offset + flat
        for i, s in zip(indices[1:n], self.strides[1:n]):
            flat = flat + (i if s == literal_one else s * i)

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
        self.indices = tuple(as_cexpr(i) for i in indices)

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

    def ce_format(self):
        s = self.array.ce_format()
        for index in self.indices:
            s += "[" + index.ce_format() + "]"
        return s

class Conditional(CExprOperator):
    __slots__ = ("condition", "true", "false")
    precedence = PRECEDENCE.CONDITIONAL

    def __init__(self, condition, true, false):
        self.condition = as_cexpr(condition)
        self.true = as_cexpr(true)
        self.false = as_cexpr(false)

    def ce_format(self):
        # Format children
        c = self.condition.ce_format()
        t = self.true.ce_format()
        f = self.false.ce_format()

        # Apply parentheses
        if self.condition.precedence >= self.precedence:
            c = '(' + c + ')'
        if self.true.precedence >= self.precedence:
            t = '(' + t + ')'
        if self.false.precedence >= self.precedence:
            f = '(' + f + ')'

        # Return combined string
        return c + " ? " + t + " : " + f

class Call(CExprOperator):
    __slots__ = ("function", "arguments")
    precedence = PRECEDENCE.CALL
    sideeffect = True

    def __init__(self, function, arguments=None):
        # Note: This will wrap a str as a Symbol
        self.function = as_cexpr(function)
        # Accept None, single, or multple arguments; literals or CExprs
        if arguments is None:
            arguments = ()
        elif not isinstance(arguments, (tuple, list)):
            arguments = (arguments,)
        self.arguments = [as_cexpr(arg) for arg in arguments]

    def ce_format(self):
        args = ", ".join(arg.ce_format() for arg in self.arguments)
        return self.function.ce_format() + "(" + args + ")"


############## Convertion function to expression nodes

number_types = (int, float, complex, numpy.number)

def as_cexpr(node):
    """Typechecks and wraps an object as a valid CExpr.

    Accepts CExpr nodes, treats int and float as literals, and treats a string as a symbol.
    """
    global number_types
    if isinstance(node, CExpr):
        return node
    elif isinstance(node, (int, numpy.integer)):
        return LiteralInt(node)
    elif isinstance(node, (float, numpy.floating)):
        return LiteralFloat(node)
    elif isinstance(node, str):
        # Treat string as a symbol
        # TODO: Using LiteralString or VerbatimExpr would be other options, is this too ambiguous?
        return Symbol(node)
    else:
        raise RuntimeError("Unexpected CExpr type %s:\n%s" % (type(node), str(node)))

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
    __slots__ = ("children",)

    def cs_format(self):
        "Return S: string | list(S) | Indented(S)."
        raise NotImplementedError("Missing implementation of cs_format() in CStatement.")

    def __str__(self):
        try:
            s = self.cs_format()
        except:
            if CNode.debug:
                print("Error in CStatement string formatting. Inspect self.")
                import IPython; IPython.embed()
            raise
        return format_indented_lines(s)


############## Statements

class VerbatimStatement(CStatement):
    "Wraps a source code string to be pasted verbatim into the source code."
    __slots__ = ("codestring",)
    def __init__(self, codestring):
        assert isinstance(codestring, str)
        self.codestring = codestring

    def cs_format(self):
        return self.codestring

class Statement(CStatement):
    "Make an expression into a statement."
    __slots__ = ("expr",)
    def __init__(self, expr):
        self.expr = as_cexpr(expr)

    def cs_format(self):
        return self.expr.ce_format() + ";"

class StatementList(CStatement):
    "A simple sequence of statements. No new scopes are introduced."
    __slots__ = ("statements",)
    def __init__(self, statements):
        self.statements = [as_cstatement(st) for st in statements]

    def cs_format(self):
        return [st.cs_format() for st in self.statements]


############## Simple statements

class Using(CStatement):
    __slots__ = ("name",)
    def __init__(self, name):
        assert isinstance(name, str)
        self.name = name

    def cs_format(self):
        return "using " + self.name + ";"

class Break(CStatement):
    __slots__ = ()
    def cs_format(self):
        return "break;"

class Continue(CStatement):
    __slots__ = ()
    def cs_format(self):
        return "continue;"

class Return(CStatement):
    __slots__ = ("value",)
    def __init__(self, value):
        self.value = as_cexpr(value)

    def cs_format(self):
        return "return " + self.value.ce_format() + ";"

class Case(CStatement):
    __slots__ = ("value",)
    def __init__(self, value):
        # NB! This is too permissive and will allow invalid case arguments.
        self.value = as_cexpr(value)

    def cs_format(self):
        return "case " + self.value.ce_format() + ":"

class Default(CStatement):
    __slots__ = ()
    def cs_format(self):
        return "default:"

class Throw(CStatement):
    __slots__ = ("exception", "message")
    def __init__(self, exception, message):
        assert isinstance(exception, str)
        assert isinstance(message, str)
        self.exception = exception
        self.message = message

    def cs_format(self):
        assert '"' not in self.message
        return "throw " + self.exception + '("' + self.message + '");'

class Comment(CStatement):
    "Line comment(s) used for annotating the generated code with human readable remarks."
    __slots__ = ("comment",)
    def __init__(self, comment):
        assert isinstance(comment, str)
        self.comment = comment

    def cs_format(self):
        lines = self.comment.strip().split("\n")
        return ["// " + line.strip() for line in lines]

class Pragma(CStatement): # TODO: Improve on this with a use case later
    "Pragma comments used for compiler-specific annotations."
    __slots__ = ("comment",)
    def __init__(self, comment):
        assert isinstance(comment, str)
        self.comment = comment

    def cs_format(self):
        assert "\n" not in self.comment
        return "#pragma " + self.comment


############## Type and variable declarations

class VariableDecl(CStatement):
    "Declare a variable, optionally define initial value."
    __slots__ = ("typename", "symbol", "value")
    def __init__(self, typename, symbol, value=None):

        # No type system yet, just using strings
        assert isinstance(typename, str)
        self.typename = typename

        # Allow Symbol or just a string
        self.symbol = as_symbol(symbol)

        if value is not None:
            value = as_cexpr(value)
        self.value = value

    def cs_format(self):
        code = self.typename + " " + self.symbol.name
        if self.value is not None:
            code += " = " + self.value.ce_format()
        return code + ";"


def leftover(size, padlen):
    "Return minimum integer to add to size to make it divisible by padlen."
    return (padlen - (size % padlen)) % padlen


def build_1d_initializer_list(values, formatter, padlen=0):
    '''Return a list containing a single line formatted like "{ 0.0, 1.0, 2.0 }"'''
    tokens = ["{ "]
    if numpy.product(values.shape) > 0:
        sep = ", "
        fvalues = [formatter(v) for v in values]
        for v in fvalues[:-1]:
            tokens.append(v)
            tokens.append(sep)
        tokens.append(fvalues[-1])
        if padlen:
            # Add padding
            zero = formatter(0)
            for i in range(leftover(len(values), padlen)):
                tokens.append(sep)
                tokens.append(zero)
    tokens += " }"
    return "".join(tokens)


def build_initializer_lists(values, sizes, level, formatter, padlen=0):
    """Return a list of lines with initializer lists for a multidimensional array.

    Example output::

        { { 0.0, 0.1 },
          { 1.0, 1.1 } }
    """
    values = numpy.asarray(values)
    assert numpy.product(values.shape) == numpy.product(sizes)
    assert len(sizes) > 0
    assert len(values.shape) > 0
    assert len(sizes) == len(values.shape)
    assert numpy.all(values.shape == sizes)

    r = len(sizes)
    assert r > 0
    if r == 1:
        return [build_1d_initializer_list(values, formatter, padlen=padlen)]
    else:
        # Render all sublists
        parts = []
        for val in values:
            sublist = build_initializer_lists(val, sizes[1:], level+1, formatter, padlen=padlen)
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

def _is_zero(values):
    if isinstance(values, (int, float, LiteralFloat, LiteralInt)):
        return float(values) == 0.0
    else:
        return numpy.count_nonzero(values) == 0

class ArrayDecl(CStatement):
    """A declaration or definition of an array.

    Note that just setting values=0 is sufficient
    to initialize the entire array to zero.

    Otherwise use nested lists of lists to represent
    multidimensional array values to initialize to.
    """
    __slots__ = ("typename", "symbol", "sizes", "alignas", "padlen", "values")
    def __init__(self, typename, symbol, sizes, values=None, alignas=None, padlen=0):
        assert isinstance(typename, str)
        self.typename = typename

        self.symbol = as_symbol(symbol)

        if isinstance(sizes, int):
            sizes = (sizes,)
        self.sizes = tuple(sizes)

        # NB! No type checking, assuming nested lists of literal values. Not applying as_cexpr.
        self.values = values

        self.alignas = alignas
        self.padlen = padlen

    def __getitem__(self, indices):
        """Allow using array declaration object as the array when indexed.

        A = ArrayDecl("int", "A", (2,3))
        code = [A, Assign(A[0,0], 1.0)]
        """
        return ArrayAccess(self, indices)

    def cs_format(self):
        # Pad innermost array dimension
        sizes = list(self.sizes)
        if self.padlen:
            sizes[-1] += leftover(sizes[-1], self.padlen)

        # C style
        brackets = ''.join("[%d]" % n for n in sizes)
        decl = self.typename + " " + self.symbol.name + brackets

        # C++11 style with std::array # TODO: Enable this, needs #include <array>
        #typename = self.typename
        #for dim in reversed(sizes):
        #    typename = "std::array<%s, %s>" % (typename, dim)
        #decl = "%s %s" % (typename, self.symbol.name)

        # TODO: support aligning inner dimensions with padding
        # C++11 style alignas prefix
        if self.alignas:
            align = "alignas(%d)" % int(self.alignas)
            decl = align + " " + decl

        if self.values is None:
            # Undefined initial values
            return decl + ";"
        elif _is_zero(self.values):
            # Zero initial values
            return decl + " = {};"
        else:
            # Construct initializer lists for arbitrary multidimensional array values
            initializer_lists = build_initializer_lists(self.values, self.sizes, 0,
                                                        format_value, padlen=self.padlen)
            if len(initializer_lists) == 1:
                return decl + " = " + initializer_lists[0] + ";"
            else:
                initializer_lists[-1] += ";" # Close statement on final line
                return (decl + " =", Indented(initializer_lists))


############## Scoped statements

class Scope(CStatement):
    __slots__ = ("body",)
    def __init__(self, body):
        self.body = as_cstatement(body)

    def cs_format(self):
        return ("{", Indented(self.body.cs_format()), "}")

class Namespace(CStatement):
    __slots__ = ("name", "body")
    def __init__(self, name, body):
        assert isinstance(name, str)
        self.name = name
        self.body = as_cstatement(body)

    def cs_format(self):
        return ("namespace " + self.name,
                "{", Indented(self.body.cs_format()), "}")

class If(CStatement):
    __slots__ = ("condition", "body")
    def __init__(self, condition, body):
        self.condition = as_cexpr(condition)
        self.body = as_cstatement(body)

    def cs_format(self):
        return ("if (" + self.condition.ce_format() + ")",
                "{", Indented(self.body.cs_format()), "}")

class ElseIf(CStatement):
    __slots__ = ("condition", "body")
    def __init__(self, condition, body):
        self.condition = as_cexpr(condition)
        self.body = as_cstatement(body)

    def cs_format(self):
        return ("else if (" + self.condition.ce_format() + ")",
                "{", Indented(self.body.cs_format()), "}")

class Else(CStatement):
    __slots__ = ("body",)
    def __init__(self, body):
        self.body = as_cstatement(body)

    def cs_format(self):
        return ("else",
                "{", Indented(self.body.cs_format()), "}")

class While(CStatement):
    __slots__ = ("condition", "body")
    def __init__(self, condition, body):
        self.condition = as_cexpr(condition)
        self.body = as_cstatement(body)

    def cs_format(self):
        return ("while (" + self.condition.ce_format() + ")",
                "{", Indented(self.body.cs_format()), "}")

class Do(CStatement):
    __slots__ = ("condition", "body")
    def __init__(self, condition, body):
        self.condition = as_cexpr(condition)
        self.body = as_cstatement(body)

    def cs_format(self):
        return ("do", "{", Indented(self.body.cs_format()),
                "} while (" + self.condition.ce_format() + ");")

class For(CStatement):
    __slots__ = ("init", "check", "update", "body")
    def __init__(self, init, check, update, body):
        self.init = as_cstatement(init)
        self.check = as_cexpr(check)
        self.update = as_cexpr(update)
        self.body = as_cstatement(body)

    def cs_format(self):
        # The C model here is a bit crude and this causes trouble
        # in the init statement/expression here:
        init = self.init.cs_format()
        assert isinstance(init, str)
        assert init.rstrip().endswith(";")

        check = self.check.ce_format()
        update = self.update.ce_format()
        body = self.body.cs_format()
        return ("for (" + init + " " + check + "; " + update + ")",
                "{", Indented(body), "}")

class Switch(CStatement):
    __slots__ = ("arg", "cases", "default", "autobreak", "autoscope")
    def __init__(self, arg, cases, default=None, autobreak=True, autoscope=True):
        self.arg = as_cexpr(arg)
        self.cases = [(as_cexpr(value), as_cstatement(body)) for value, body in cases]
        if default is not None:
            default = as_cstatement(default)
        self.default = default
        assert autobreak in (True, False)
        assert autoscope in (True, False)
        self.autobreak = autobreak
        self.autoscope = autoscope

    def cs_format(self):
        cases = []
        for case in self.cases:
            caseheader = "case " + case[0].ce_format() + ":"
            casebody = case[1].cs_format()
            if self.autoscope:
                casebody = ("{", Indented(casebody), "}")
            if self.autobreak:
                casebody = (casebody, "break;")
            cases.extend([caseheader, Indented(casebody)])

        if self.default is not None:
            caseheader = "default:"
            casebody = self.default.cs_format()
            if self.autoscope:
                casebody = ("{", Indented(casebody), "}")
            cases.extend([caseheader, Indented(casebody)])

        return ("switch (" + self.arg.ce_format() + ")",
                "{", cases, "}")


class ForRange(CStatement):
    "Slightly higher-level for loop assuming incrementing an index over a range."
    __slots__ = ("index", "begin", "end", "body", "index_type")
    def __init__(self, index, begin, end, body):
        self.index = as_cexpr(index)
        self.begin = as_cexpr(begin)
        self.end = as_cexpr(end)
        self.body = as_cstatement(body)

        # Could be configured if needed but not sure how we're
        # going to handle type information right now:
        self.index_type = "int"

    def cs_format(self):
        indextype = self.index_type
        index = self.index.ce_format()
        begin = self.begin.ce_format()
        end = self.end.ce_format()

        init = indextype + " " + index + " = " + begin
        check = index + " < " + end
        update = "++" + index

        return ("for (" + init + "; " + check + "; " + update + ")",
                "{", Indented(self.body.cs_format()), "}")


############## Convertion function to statement nodes

def as_cstatement(node):
    "Perform type checking on node and wrap in a suitable statement type if necessary."
    if isinstance(node, CStatement):
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
        return StatementList(node)
    elif isinstance(node, str):
        # Backdoor for flexibility in code generation to allow verbatim pasted statements
        return VerbatimStatement(node)
    else:
        raise RuntimeError("Unexpected CStatement type %s:\n%s" % (type(node), str(node)))
