# Copyright (C) 2023 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later


math_table = {
    "double": {
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
    },
    "float": {
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
    },
    "long double": {
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
    "double _Complex": {
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
    },
    "float _Complex": {
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
    },
}


def build_initializer_lists(values):
    arr = "{"
    if len(values.shape) == 1:
        arr += ", ".join(str(v) for v in values)
    elif len(values.shape) > 1:
        arr += ", \n".join(build_initializer_lists(v) for v in values)
    arr += "}"
    return arr


class CFormatter(object):
    def __init__(self, scalar) -> None:
        self.scalar_type = scalar

    def format_statement_list(self, slist) -> str:
        output = ""
        for s in slist.statements:
            output += self.c_format(s)
        return output

    def format_comment(self, c) -> str:
        return "// " + c.comment + "\n"

    def format_array_decl(self, arr) -> str:
        symbol = self.c_format(arr.symbol)
        dims = "".join([f"[{i}]" for i in arr.sizes])
        if arr.values is None:
            vals = "{}"
        else:
            vals = build_initializer_lists(arr.values)
        return f"{arr.typename} {symbol}{dims} = {vals};\n"

    def format_array_access(self, arr) -> str:
        name = self.c_format(arr.array)
        indices = f"[{']['.join(self.c_format(i) for i in arr.indices)}]"
        return f"{name}{indices}"

    def format_variable_decl(self, v) -> str:
        val = self.c_format(v.value)
        symbol = self.c_format(v.symbol)
        return f"{v.typename} {symbol} = {val};\n"

    def format_nary_op(self, oper) -> str:
        # Format children
        args = [self.c_format(arg) for arg in oper.args]

        # Apply parentheses
        for i in range(len(args)):
            if oper.args[i].precedence >= oper.precedence:
                args[i] = "(" + args[i] + ")"

        # Return combined string
        return f" {oper.op} ".join(args)

    def format_binary_op(self, oper) -> str:
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

    def format_literal_float(self, val) -> str:
        return f"{val.value}"

    def format_literal_int(self, val) -> str:
        return f"{val.value}"

    def format_for_range(self, r) -> str:
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
        return self.c_format(s.expr)

    def format_assign(self, expr) -> str:
        rhs = self.c_format(expr.rhs)
        lhs = self.c_format(expr.lhs)
        return f"{lhs} {expr.op} {rhs};\n"

    def format_conditional(self, s) -> str:
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
        return f"{s.name}"

    def format_math_function(self, c) -> str:
        # A few translations to get tests working - need to do properly.
        # Depending on dtype...
        math_table = {"power": "pow", "abs": "fabs"}
        if "_Complex" in self.scalar_type:
            math_table = {
                "power": "cpow",
                "real": "creal",
                "imag": "cimag",
                "abs": "cabs",
            }
        func = math_table.get(c.function, c.function)
        args = ", ".join(self.c_format(arg) for arg in c.args)
        return f"{func}({args})"

    c_impl = {
        "StatementList": format_statement_list,
        "Comment": format_comment,
        "ArrayDecl": format_array_decl,
        "ArrayAccess": format_array_access,
        "VariableDecl": format_variable_decl,
        "ForRange": format_for_range,
        "Statement": format_statement,
        "Assign": format_assign,
        "AssignAdd": format_assign,
        "Product": format_nary_op,
        "Sum": format_nary_op,
        "Add": format_binary_op,
        "Mul": format_binary_op,
        "Div": format_binary_op,
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
        name = s.__class__.__name__
        try:
            return self.c_impl[name](self, s)
        except KeyError:
            raise RuntimeError("Unknown statement: ", name)
