# Copyright (C) 2023 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import ffcx.codegeneration.lnodes as L

math_table = {
    "sqrt": "std::sqrt",
    "abs": "std::abs",
    "cos": "std::cos",
    "sin": "std::sin",
    "tan": "std::tan",
    "acos": "std::acos",
    "asin": "std::asin",
    "atan": "std::atan",
    "cosh": "std::cosh",
    "sinh": "std::sinh",
    "tanh": "std::tanh",
    "acosh": "std::acosh",
    "asinh": "std::asinh",
    "atanh": "std::atanh",
    "power": "std::pow",
    "exp": "std::exp",
    "ln": "std::log",
    "erf": "std::erf",
    "atan_2": "std::atan2",
    "min_value": "std::fmin",
    "max_value": "std::fmax",
    "bessel_y": "std::cyl_bessel_i",
    "bessel_j": "std::cyl_bessel_j",
    "conj": "std::conj",
    "real": "std::real",
    "imag": "std::imag"}


def build_initializer_lists(values):
    arr = "{"
    if len(values.shape) == 1:
        return "{" + ", ".join(str(v) for v in values) + "}"
    elif len(values.shape) > 1:
        arr += ",\n".join(build_initializer_lists(v) for v in values)
    arr += "}"
    return arr


class CppFormatter(object):
    def __init__(self, scalar) -> None:
        self.scalar_type = "T"
        self.real_type = "U"

    def format_statement_list(self, slist) -> str:
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
        return "// " + c.comment + "\n"

    def format_array_decl(self, arr) -> str:
        dtype = arr.symbol.dtype
        assert dtype is not None

        if dtype == L.DataType.SCALAR:
            typename = self.scalar_type
        elif dtype == L.DataType.REAL:
            typename = self.real_type
        elif dtype == L.DataType.INT:
            typename = "int"
        else:
            raise ValueError("Invalid datatype")

        symbol = self.c_format(arr.symbol)
        dims = "".join([f"[{i}]" for i in arr.sizes])
        if arr.values is None:
            assert arr.const is False
            return f"{typename} {symbol}{dims};\n"

        vals = build_initializer_lists(arr.values)
        cstr = "static const " if arr.const else ""
        return f"{cstr}{typename} {symbol}{dims} = {vals};\n"

    def format_array_access(self, arr) -> str:
        name = self.c_format(arr.array)
        indices = f"[{']['.join(self.c_format(i) for i in arr.indices)}]"
        return f"{name}{indices}"

    def format_multi_index(self, index) -> str:
        return self.c_format(index.global_index)

    def format_variable_decl(self, v) -> str:
        val = self.c_format(v.value)
        symbol = self.c_format(v.symbol)
        assert v.symbol.dtype
        if v.symbol.dtype == L.DataType.SCALAR:
            typename = self.scalar_type
        elif v.symbol.dtype == L.DataType.REAL:
            typename = self.real_type
        return f"{typename} {symbol} = {val};\n"

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

    def format_neg(self, val) -> str:
        arg = self.c_format(val.arg)
        return f"-{arg}"

    def format_not(self, val) -> str:
        arg = self.c_format(val.arg)
        return f"{val.op}({arg})"

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
        # Get a function from the table, if available, else just use bare name
        func = math_table.get(c.function, c.function)
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
        "Neg": format_neg,
        "Sum": format_nary_op,
        "Add": format_binary_op,
        "Sub": format_binary_op,
        "Mul": format_binary_op,
        "Div": format_binary_op,
        "Not": format_not,
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
