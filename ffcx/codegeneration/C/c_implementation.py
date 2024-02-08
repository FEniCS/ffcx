# Copyright (C) 2023 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""C implementation."""

import warnings

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


class CFormatter:
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

    def format_statement_list(self, slist) -> str:
        """Format a statement list."""
        return "".join(self.c_format(s) for s in slist.statements)

    def format_section(self, section) -> str:
        """Format a section."""
        # add new line before section
        comments = "// ------------------------ \n"
        comments += "// Section: " + section.name + "\n"
        comments += "// Inputs: " + ", ".join(w.name for w in section.input) + "\n"
        comments += "// Outputs: " + ", ".join(w.name for w in section.output) + "\n"
        declarations = "".join(self.c_format(s) for s in section.declarations)

        body = ""
        if len(section.statements) > 0:
            declarations += "{\n  "
            body = "".join(self.c_format(s) for s in section.statements)
            body = body.replace("\n", "\n  ")
            body = body[:-2] + "}\n"

        body += "// ------------------------ \n"
        return comments + declarations + body

    def format_comment(self, c) -> str:
        """Format a comment."""
        return "// " + c.comment + "\n"

    def format_array_decl(self, arr) -> str:
        """Format an array declaration."""
        dtype = arr.symbol.dtype
        typename = self._dtype_to_name(dtype)

        symbol = self.c_format(arr.symbol)
        dims = "".join([f"[{i}]" for i in arr.sizes])
        if arr.values is None:
            assert arr.const is False
            return f"{typename} {symbol}{dims};\n"

        vals = self._build_initializer_lists(arr.values)
        cstr = "static const " if arr.const else ""
        return f"{cstr}{typename} {symbol}{dims} = {vals};\n"

    def format_array_access(self, arr) -> str:
        """Format an array access."""
        name = self.c_format(arr.array)
        indices = f"[{']['.join(self.c_format(i) for i in arr.indices)}]"
        return f"{name}{indices}"

    def format_variable_decl(self, v) -> str:
        """Format a variable declaration."""
        val = self.c_format(v.value)
        symbol = self.c_format(v.symbol)
        typename = self._dtype_to_name(v.symbol.dtype)
        return f"{typename} {symbol} = {val};\n"

    def format_nary_op(self, oper) -> str:
        """Format an n-ary operation."""
        # Format children
        args = [self.c_format(arg) for arg in oper.args]

        # Apply parentheses
        for i in range(len(args)):
            if oper.args[i].precedence >= oper.precedence:
                args[i] = "(" + args[i] + ")"

        # Return combined string
        return f" {oper.op} ".join(args)

    def format_binary_op(self, oper) -> str:
        """Format a binary operation."""
        # Format children
        lhs = self.c_format(oper.lhs)
        rhs = self.c_format(oper.rhs)

        # Apply parentheses
        if oper.lhs.precedence >= oper.precedence:
            lhs = f"({lhs})"
        if oper.rhs.precedence >= oper.precedence:
            rhs = f"({rhs})"

        # Return combined string
        return f"{lhs} {oper.op} {rhs}"

    def format_unary_op(self, oper) -> str:
        """Format a unary operation."""
        arg = self.c_format(oper.arg)
        if oper.arg.precedence >= oper.precedence:
            return f"{oper.op}({arg})"
        return f"{oper.op}{arg}"

    def format_literal_float(self, val) -> str:
        """Format a literal float."""
        value = self._format_number(val.value)
        return f"{value}"

    def format_literal_int(self, val) -> str:
        """Format a literal int."""
        return f"{val.value}"

    def format_for_range(self, r) -> str:
        """Format a for loop over a range."""
        begin = self.c_format(r.begin)
        end = self.c_format(r.end)
        index = self.c_format(r.index)
        output = f"for (int {index} = {begin}; {index} < {end}; ++{index})\n"
        output += "{\n"
        body = self.c_format(r.body)
        for line in body.split("\n"):
            if len(line) > 0:
                output += f"  {line}\n"
        output += "}\n"
        return output

    def format_statement(self, s) -> str:
        """Format a statement."""
        return self.c_format(s.expr)

    def format_assign(self, expr) -> str:
        """Format an assignment."""
        rhs = self.c_format(expr.rhs)
        lhs = self.c_format(expr.lhs)
        return f"{lhs} {expr.op} {rhs};\n"

    def format_conditional(self, s) -> str:
        """Format a conditional."""
        # Format children
        c = self.c_format(s.condition)
        t = self.c_format(s.true)
        f = self.c_format(s.false)

        # Apply parentheses
        if s.condition.precedence >= s.precedence:
            c = "(" + c + ")"
        if s.true.precedence >= s.precedence:
            t = "(" + t + ")"
        if s.false.precedence >= s.precedence:
            f = "(" + f + ")"

        # Return combined string
        return c + " ? " + t + " : " + f

    def format_symbol(self, s) -> str:
        """Format a symbol."""
        return f"{s.name}"

    def format_multi_index(self, mi) -> str:
        """Format a multi-index."""
        return self.c_format(mi.global_index)

    def format_math_function(self, c) -> str:
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
        args = ", ".join(self.c_format(arg) for arg in c.args)
        return f"{func}({args})"

    c_impl = {
        "Section": format_section,
        "StatementList": format_statement_list,
        "Comment": format_comment,
        "ArrayDecl": format_array_decl,
        "ArrayAccess": format_array_access,
        "MultiIndex": format_multi_index,
        "VariableDecl": format_variable_decl,
        "ForRange": format_for_range,
        "Statement": format_statement,
        "Assign": format_assign,
        "AssignAdd": format_assign,
        "Product": format_nary_op,
        "Neg": format_unary_op,
        "Sum": format_nary_op,
        "Add": format_binary_op,
        "Sub": format_binary_op,
        "Mul": format_binary_op,
        "Div": format_binary_op,
        "Not": format_unary_op,
        "LiteralFloat": format_literal_float,
        "LiteralInt": format_literal_int,
        "Symbol": format_symbol,
        "Conditional": format_conditional,
        "MathFunction": format_math_function,
        "And": format_binary_op,
        "Or": format_binary_op,
        "NE": format_binary_op,
        "EQ": format_binary_op,
        "GE": format_binary_op,
        "LE": format_binary_op,
        "GT": format_binary_op,
        "LT": format_binary_op,
    }

    def c_format(self, s) -> str:
        """Format as C."""
        name = s.__class__.__name__
        try:
            return self.c_impl[name](self, s)
        except KeyError:
            raise RuntimeError("Unknown statement: ", name)
