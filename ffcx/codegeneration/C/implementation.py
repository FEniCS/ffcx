# Copyright (C) 2023-2025 Chris Richardson and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""C implementation."""

import warnings
from functools import singledispatchmethod

import numpy as np
import numpy.typing as npt

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype

math_table = {
    "float64": {
        "sqrt": "sqrt",
        "abs": "fabs",
        "cos": "cos",
        "sin": "sin",
        "tan": "tan",
        "acos": "acos",
        "asin": "asin",
        "atan": "atan",
        "cosh": "cosh",
        "sinh": "sinh",
        "tanh": "tanh",
        "acosh": "acosh",
        "asinh": "asinh",
        "atanh": "atanh",
        "power": "pow",
        "exp": "exp",
        "ln": "log",
        "erf": "erf",
        "atan_2": "atan2",
        "min_value": "fmin",
        "max_value": "fmax",
        "bessel_y": "yn",
        "bessel_j": "jn",
    },
    "float32": {
        "sqrt": "sqrtf",
        "abs": "fabsf",
        "cos": "cosf",
        "sin": "sinf",
        "tan": "tanf",
        "acos": "acosf",
        "asin": "asinf",
        "atan": "atanf",
        "cosh": "coshf",
        "sinh": "sinhf",
        "tanh": "tanhf",
        "acosh": "acoshf",
        "asinh": "asinhf",
        "atanh": "atanhf",
        "power": "powf",
        "exp": "expf",
        "ln": "logf",
        "erf": "erff",
        "atan_2": "atan2f",
        "min_value": "fminf",
        "max_value": "fmaxf",
        "bessel_y": "yn",
        "bessel_j": "jn",
    },
    "longdouble": {
        "sqrt": "sqrtl",
        "abs": "fabsl",
        "cos": "cosl",
        "sin": "sinl",
        "tan": "tanl",
        "acos": "acosl",
        "asin": "asinl",
        "atan": "atanl",
        "cosh": "coshl",
        "sinh": "sinhl",
        "tanh": "tanhl",
        "acosh": "acoshl",
        "asinh": "asinhl",
        "atanh": "atanhl",
        "power": "powl",
        "exp": "expl",
        "ln": "logl",
        "erf": "erfl",
        "atan_2": "atan2l",
        "min_value": "fminl",
        "max_value": "fmaxl",
    },
    "complex128": {
        "sqrt": "csqrt",
        "abs": "cabs",
        "cos": "ccos",
        "sin": "csin",
        "tan": "ctan",
        "acos": "cacos",
        "asin": "casin",
        "atan": "catan",
        "cosh": "ccosh",
        "sinh": "csinh",
        "tanh": "ctanh",
        "acosh": "cacosh",
        "asinh": "casinh",
        "atanh": "catanh",
        "power": "cpow",
        "exp": "cexp",
        "ln": "clog",
        "real": "creal",
        "imag": "cimag",
        "conj": "conj",
        "max_value": "fmax",
        "min_value": "fmin",
        "bessel_y": "yn",
        "bessel_j": "jn",
    },
    "complex64": {
        "sqrt": "csqrtf",
        "abs": "cabsf",
        "cos": "ccosf",
        "sin": "csinf",
        "tan": "ctanf",
        "acos": "cacosf",
        "asin": "casinf",
        "atan": "catanf",
        "cosh": "ccoshf",
        "sinh": "csinhf",
        "tanh": "ctanhf",
        "acosh": "cacoshf",
        "asinh": "casinhf",
        "atanh": "catanhf",
        "power": "cpowf",
        "exp": "cexpf",
        "ln": "clogf",
        "real": "crealf",
        "imag": "cimagf",
        "conj": "conjf",
        "max_value": "fmaxf",
        "min_value": "fminf",
        "bessel_y": "yn",
        "bessel_j": "jn",
    },
}


class Formatter:
    """C formatter."""

    scalar_type: np.dtype
    real_type: np.dtype

    def __init__(self, dtype: npt.DTypeLike) -> None:
        """Initialise."""
        self.scalar_type = np.dtype(dtype)
        self.real_type = dtype_to_scalar_dtype(dtype)

    def _dtype_to_name(self, dtype) -> str:
        """Convert dtype to C name."""
        if dtype == L.DataType.SCALAR:
            return dtype_to_c_type(self.scalar_type)
        if dtype == L.DataType.REAL:
            return dtype_to_c_type(self.real_type)
        if dtype == L.DataType.INT:
            return "int"
        if dtype == L.DataType.BOOL:
            return "bool"
        raise ValueError(f"Invalid dtype: {dtype}")

    def _format_number(self, x):
        """Format a number."""
        # Use 16sf for precision (good for float64 or less)
        if isinstance(x, complex):
            return f"({x.real:.16}+I*{x.imag:.16})"
        elif isinstance(x, float):
            return f"{x:.16}"
        return str(x)

    def _build_initializer_lists(self, values):
        """Build initializer lists."""
        arr = "{"
        if len(values.shape) == 1:
            arr += ", ".join(self._format_number(v) for v in values)
        elif len(values.shape) > 1:
            arr += ",\n  ".join(self._build_initializer_lists(v) for v in values)
        arr += "}"
        return arr

    @singledispatchmethod
    def __call__(self, obj: L.LNode) -> str:
        """Format an L Node."""
        raise NotImplementedError(f"Can not format object to type {type(obj)}")

    @__call__.register
    def _(self, slist: L.StatementList) -> str:
        """Format a statement list."""
        return "".join(self(s) for s in slist.statements)

    @__call__.register
    def _(self, section: L.Section) -> str:
        """Format a section."""
        # add new line before section
        comments = (
            f"// ------------------------ \n"
            f"// Section: {section.name}\n"
            f"// Inputs: {', '.join(w.name for w in section.input)}\n"
            f"// Outputs: {', '.join(w.name for w in section.output)}\n"
        )
        declarations = "".join(self(s) for s in section.declarations)

        body = ""
        if len(section.statements) > 0:
            declarations += "{\n  "
            body = "".join(self(s) for s in section.statements)
            body = body.replace("\n", "\n  ")
            body = body[:-2] + "}\n"

        body += "// ------------------------ \n"
        return str(comments + declarations + body)

    @__call__.register
    def _(self, c: L.Comment) -> str:
        """Format a comment."""
        return "// " + c.comment + "\n"

    @__call__.register
    def _(self, arr: L.ArrayDecl) -> str:
        """Format an array declaration."""
        dtype = arr.symbol.dtype
        typename = self._dtype_to_name(dtype)

        symbol = self(arr.symbol)
        dims = "".join([f"[{i}]" for i in arr.sizes])
        if arr.values is None:
            assert arr.const is False
            return f"{typename} {symbol}{dims};\n"

        vals = self._build_initializer_lists(arr.values)
        cstr = "static const " if arr.const else ""
        return f"{cstr}{typename} {symbol}{dims} = {vals};\n"

    @__call__.register
    def _(self, arr: L.ArrayAccess) -> str:
        """Format an array access."""
        name = self(arr.array)
        indices = f"[{']['.join(self(i) for i in arr.indices)}]"
        return f"{name}{indices}"

    @__call__.register
    def _(self, v: L.VariableDecl) -> str:
        """Format a variable declaration."""
        val = self(v.value)
        symbol = self(v.symbol)
        typename = self._dtype_to_name(v.symbol.dtype)
        return f"{typename} {symbol} = {val};\n"

    @__call__.register
    def _(self, oper: L.NaryOp) -> str:
        """Format an n-ary operation."""
        # Format children
        args = [self(arg) for arg in oper.args]

        # Apply parentheses
        for i in range(len(args)):
            if oper.args[i].precedence >= oper.precedence:
                args[i] = "(" + args[i] + ")"

        # Return combined string
        return f" {oper.op} ".join(args)

    @__call__.register
    def _(self, oper: L.BinOp) -> str:
        """Format a binary operation."""
        # Format children
        lhs = self(oper.lhs)
        rhs = self(oper.rhs)

        # Apply parentheses
        if oper.lhs.precedence >= oper.precedence:
            lhs = f"({lhs})"
        if oper.rhs.precedence >= oper.precedence:
            rhs = f"({rhs})"

        # Return combined string
        return f"{lhs} {oper.op} {rhs}"

    @__call__.register(L.Neg)
    @__call__.register(L.Not)
    def _(self, oper) -> str:
        """Format a unary operation."""
        arg = self(oper.arg)
        if oper.arg.precedence >= oper.precedence:
            return f"{oper.op}({arg})"
        return f"{oper.op}{arg}"

    @__call__.register
    def _(self, val: L.LiteralFloat) -> str:
        """Format a literal float."""
        value = self._format_number(val.value)
        return f"{value}"

    @__call__.register
    def _(self, val: L.LiteralInt) -> str:
        """Format a literal int."""
        return f"{val.value}"

    @__call__.register
    def _(self, r: L.ForRange) -> str:
        """Format a for loop over a range."""
        begin = self(r.begin)
        end = self(r.end)
        index = self(r.index)
        output = f"for (int {index} = {begin}; {index} < {end}; ++{index})\n"
        output += "{\n"
        body = self(r.body)
        for line in body.split("\n"):
            if len(line) > 0:
                output += f"  {line}\n"
        output += "}\n"
        return output

    @__call__.register
    def _(self, s: L.Statement) -> str:
        """Format a statement."""
        return self(s.expr)

    @__call__.register(L.Assign)
    @__call__.register(L.AssignAdd)
    def _(self, expr: L.Assign | L.AssignAdd) -> str:
        """Format an assignment."""
        rhs = self(expr.rhs)
        lhs = self(expr.lhs)
        return f"{lhs} {expr.op} {rhs};\n"

    @__call__.register
    def _(self, s: L.Conditional) -> str:
        """Format a conditional."""
        # Format children
        c = self(s.condition)
        t = self(s.true)
        f = self(s.false)

        # Apply parentheses
        if s.condition.precedence >= s.precedence:
            c = "(" + c + ")"
        if s.true.precedence >= s.precedence:
            t = "(" + t + ")"
        if s.false.precedence >= s.precedence:
            f = "(" + f + ")"

        # Return combined string
        return c + " ? " + t + " : " + f

    @__call__.register
    def _(self, s: L.Symbol) -> str:
        """Format a symbol."""
        return f"{s.name}"

    @__call__.register
    def _(self, mi: L.MultiIndex) -> str:
        """Format a multi-index."""
        return self(mi.global_index)

    @__call__.register
    def _(self, c: L.MathFunction) -> str:
        """Format a mathematical function."""
        # Get a table of functions for this type, if available
        arg_type = self.scalar_type
        if hasattr(c.args[0], "dtype"):
            if c.args[0].dtype == L.DataType.REAL:
                arg_type = self.real_type
        else:
            warnings.warn(f"Syntax item without dtype {c.args[0]}")

        dtype_math_table = math_table[arg_type.name]

        # Get a function from the table, if available, else just use bare name
        func = dtype_math_table.get(c.function, c.function)
        args = ", ".join(self(arg) for arg in c.args)
        return f"{func}({args})"
