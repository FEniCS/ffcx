# Copyright (C) 2025 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Numba implementation for output."""

import ffcx.codegeneration.lnodes as L


def build_initializer_lists(values):
    """Build list of values."""
    arr = "["
    if len(values.shape) == 1:
        return "[" + ", ".join(str(v) for v in values) + "]"
    elif len(values.shape) > 1:
        arr += ",\n".join(build_initializer_lists(v) for v in values)
    arr += "]"
    return arr


class Formatter:
    """Implementation for numba output backend."""

    def __init__(self, scalar) -> None:
        """Initialise."""
        self.scalar_type = scalar

    def format_section(self, section: L.Section) -> str:
        """Format a section."""
        # add new line before section
        comments = "# ------------------------ \n"
        comments += "# Section: " + section.name + "\n"
        comments += "# Inputs: " + ", ".join(w.name for w in section.input) + "\n"
        comments += "# Outputs: " + ", ".join(w.name for w in section.output) + "\n"
        declarations = "".join(self.format(s) for s in section.declarations)

        body = ""
        if len(section.statements) > 0:
            body = "".join(self.format(s) for s in section.statements)

        body += "# ------------------------ \n"
        return comments + declarations + body

    def format_statement_list(self, slist: L.StatementList) -> str:
        """Format a list of statements."""
        output = ""
        for s in slist.statements:
            output += self.format(s)
        return output

    def format_comment(self, c: L.Comment) -> str:
        """Format a comment."""
        return "# " + c.comment + "\n"

    def format_array_decl(self, arr: L.ArrayDecl) -> str:
        """Format an array declaration."""
        if arr.symbol.dtype == L.DataType.SCALAR:
            dtype = "A.dtype"
        elif arr.symbol.dtype == L.DataType.REAL:
            dtype = "coordinate_dofs.dtype"
        elif arr.symbol.dtype == L.DataType.INT:
            dtype = "np.int32"
        symbol = self.format(arr.symbol)
        if arr.values is None:
            return f"{symbol} = np.empty({arr.sizes}, dtype={dtype})\n"
        elif arr.values.size == 1:
            return f"{symbol} = np.full({arr.sizes}, {arr.values[0]}, dtype={dtype})\n"
        av = build_initializer_lists(arr.values)
        av = "np.array(" + av + f", dtype={dtype})"
        return f"{symbol} = {av}\n"

    def format_array_access(self, arr: L.ArrayAccess) -> str:
        """Format array access."""
        array = self.format(arr.array)
        idx = ", ".join(self.format(ix) for ix in arr.indices)
        return f"{array}[{idx}]"

    def format_multi_index(self, index: L.MultiIndex) -> str:
        """Format a multi-index."""
        return self.format(index.global_index)

    def format_variable_decl(self, v: L.VariableDecl) -> str:
        """Format a variable declaration."""
        sym = self.format(v.symbol)
        val = self.format(v.value)
        return f"{sym} = {val}\n"

    def format_nary_op(self, oper: L.NaryOp) -> str:
        """Format a n argument operation."""
        # Format children
        args = [self.format(arg) for arg in oper.args]

        # Apply parentheses
        for i in range(len(args)):
            if oper.args[i].precedence >= oper.precedence:
                args[i] = "(" + args[i] + ")"

        # Return combined string
        return f" {oper.op} ".join(args)

    def format_binary_op(self, oper: L.BinOp) -> str:
        """Format a binary operation."""
        # Format children
        lhs = self.format(oper.lhs)
        rhs = self.format(oper.rhs)

        # Apply parentheses
        if oper.lhs.precedence >= oper.precedence:
            lhs = f"({lhs})"
        if oper.rhs.precedence >= oper.precedence:
            rhs = f"({rhs})"

        # Return combined string
        return f"{lhs} {oper.op} {rhs}"

    def format_neg(self, val: L.Neg) -> str:
        """Format unary negation."""
        arg = self.format(val.arg)
        return f"-{arg}"

    def format_not(self, val: L.Not) -> str:
        """Format not operation."""
        arg = self.format(val.arg)
        return f"not({arg})"

    def format_andor(self, oper: L.And | L.Or) -> str:
        """Format and or or operation."""
        # Format children
        lhs = self.format(oper.lhs)
        rhs = self.format(oper.rhs)

        # Apply parentheses
        if oper.lhs.precedence >= oper.precedence:
            lhs = f"({lhs})"
        if oper.rhs.precedence >= oper.precedence:
            rhs = f"({rhs})"

        opstr = {"||": "or", "&&": "and"}[oper.op]

        # Return combined string
        return f"{lhs} {opstr} {rhs}"

    def format_literal_float(self, val: L.LiteralFloat) -> str:
        """Format a literal float."""
        return f"{val.value}"

    def format_literal_int(self, val: L.LiteralInt) -> str:
        """Format a literal int."""
        return f"{val.value}"

    def format_for_range(self, r: L.ForRange) -> str:
        """Format a loop over a range."""
        begin = self.format(r.begin)
        end = self.format(r.end)
        index = self.format(r.index)
        output = f"for {index} in range({begin}, {end}):\n"
        b = self.format(r.body).split("\n")
        for line in b:
            output += f"    {line}\n"
        return output

    def format_statement(self, s: L.Statement) -> str:
        """Format a statement."""
        return self.format(s.expr)

    def format_assign(self, expr: L.Assign) -> str:
        """Format assignment."""
        rhs = self.format(expr.rhs)
        lhs = self.format(expr.lhs)
        return f"{lhs} {expr.op} {rhs}\n"

    def format_conditional(self, s: L.Conditional) -> str:
        """Format a conditional."""
        # Format children
        c = self.format(s.condition)
        t = self.format(s.true)
        f = self.format(s.false)

        # Apply parentheses
        if s.condition.precedence >= s.precedence:
            c = "(" + c + ")"
        if s.true.precedence >= s.precedence:
            t = "(" + t + ")"
        if s.false.precedence >= s.precedence:
            f = "(" + f + ")"

        # Return combined string
        return f"({t} if {c} else {f})"

    def format_symbol(self, s: L.Symbol) -> str:
        """Format a symbol."""
        return f"{s.name}"

    def format_mathfunction(self, f: L.MathFunction) -> str:
        """Format a math function."""
        function_map = {
            "ln": "log",
            "acos": "arccos",
            "asin": "arcsin",
            "atan": "arctan",
            "atan2": "arctan2",
            "acosh": "arccosh",
            "asinh": "arcsinh",
            "atanh": "arctanh",
        }
        function = function_map.get(f.function, f.function)
        args = [self.format(arg) for arg in f.args]
        if "bessel" in function:
            return "0"
        if function == "erf":
            return f"math.erf({args[0]})"
        argstr = ", ".join(args)
        return f"np.{function}({argstr})"

    impl = {
        "StatementList": format_statement_list,
        "Comment": format_comment,
        "Section": format_section,
        "ArrayDecl": format_array_decl,
        "ArrayAccess": format_array_access,
        "MultiIndex": format_multi_index,
        "VariableDecl": format_variable_decl,
        "ForRange": format_for_range,
        "Statement": format_statement,
        "Assign": format_assign,
        "AssignAdd": format_assign,
        "Product": format_nary_op,
        "Sum": format_nary_op,
        "Add": format_binary_op,
        "Sub": format_binary_op,
        "Mul": format_binary_op,
        "Div": format_binary_op,
        "Neg": format_neg,
        "Not": format_not,
        "LiteralFloat": format_literal_float,
        "LiteralInt": format_literal_int,
        "Symbol": format_symbol,
        "Conditional": format_conditional,
        "MathFunction": format_mathfunction,
        "And": format_andor,
        "Or": format_andor,
        "NE": format_binary_op,
        "EQ": format_binary_op,
        "GE": format_binary_op,
        "LE": format_binary_op,
        "GT": format_binary_op,
        "LT": format_binary_op,
    }

    def format(self, s: L.LNode) -> str:
        """Format output."""
        name = s.__class__.__name__
        try:
            return self.impl[name](self, s)
        except KeyError:
            raise RuntimeError("Unknown statement: ", name)
