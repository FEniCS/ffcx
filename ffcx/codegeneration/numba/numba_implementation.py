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


class NumbaFormatter:
    """Implementation for numba output backend."""

    def __init__(self, scalar) -> None:
        """Initialise."""
        self.scalar_type = scalar

    def format_section(self, section):
        """Format a section."""
        # add new line before section
        comments = "# ------------------------ \n"
        comments += "# Section: " + section.name + "\n"
        comments += "# Inputs: " + ", ".join(w.name for w in section.input) + "\n"
        comments += "# Outputs: " + ", ".join(w.name for w in section.output) + "\n"
        declarations = "".join(self.c_format(s) for s in section.declarations)

        body = ""
        if len(section.statements) > 0:
            body = "".join(self.c_format(s) for s in section.statements)

        body += "# ------------------------ \n"
        return comments + declarations + body

    def format_statement_list(self, slist):
        """Format a list of statements."""
        output = ""
        for s in slist.statements:
            output += self.c_format(s)
        return output

    def format_comment(self, c):
        """Format a comment."""
        return "# " + c.comment + "\n"

    def format_array_decl(self, arr):
        """Format an array declaration."""
        if arr.symbol.dtype == L.DataType.SCALAR:
            dtype = "A.dtype"
        elif arr.symbol.dtype == L.DataType.REAL:
            dtype = "coordinate_dofs.dtype"
        elif arr.symbol.dtype == L.DataType.INT:
            dtype = "np.int32"
        symbol = self.c_format(arr.symbol)
        if arr.values is None:
            return f"{symbol} = np.empty({arr.sizes}, dtype={dtype})\n"
        av = build_initializer_lists(arr.values)
        av = "np.array(" + av + f", dtype={dtype})"
        return f"{symbol} = {av}\n"

    def format_array_access(self, arr):
        """Format array access."""
        array = self.c_format(arr.array)
        idx = ", ".join(self.c_format(ix) for ix in arr.indices)
        return f"{array}[{idx}]"

    def format_multi_index(self, index):
        """Format a multi-index."""
        return self.c_format(index.global_index)

    def format_variable_decl(self, v):
        """Format a variable declaration."""
        sym = self.c_format(v.symbol)
        val = self.c_format(v.value)
        return f"{sym} = {val}\n"

    def format_nary_op(self, oper):
        """Format a n argument operation."""
        # Format children
        args = [self.c_format(arg) for arg in oper.args]

        # Apply parentheses
        for i in range(len(args)):
            if oper.args[i].precedence >= oper.precedence:
                args[i] = "(" + args[i] + ")"

        # Return combined string
        return f" {oper.op} ".join(args)

    def format_binary_op(self, oper):
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

    def format_neg(self, val):
        """Format unary negation."""
        arg = self.c_format(val.arg)
        return f"-{arg}"

    def format_not(self, val):
        """Format not operation."""
        arg = self.c_format(val.arg)
        return f"not({arg})"

    def format_andor(self, oper):
        """Format and or or operation."""
        # Format children
        lhs = self.c_format(oper.lhs)
        rhs = self.c_format(oper.rhs)

        # Apply parentheses
        if oper.lhs.precedence >= oper.precedence:
            lhs = f"({lhs})"
        if oper.rhs.precedence >= oper.precedence:
            rhs = f"({rhs})"

        opstr = {"||": "or", "&&": "and"}[oper.op]

        # Return combined string
        return f"{lhs} {opstr} {rhs}"

    def format_literal_float(self, val):
        """Format a literal float."""
        return f"{val.value}"

    def format_literal_int(self, val):
        """Format a literal int."""
        return f"{val.value}"

    def format_for_range(self, r):
        """Format a loop over a range."""
        begin = self.c_format(r.begin)
        end = self.c_format(r.end)
        index = self.c_format(r.index)
        output = f"for {index} in range({begin}, {end}):\n"
        b = self.c_format(r.body).split("\n")
        for line in b:
            output += f"    {line}\n"
        return output

    def format_statement(self, s):
        """Format a statement."""
        return self.c_format(s.expr)

    def format_assign(self, expr):
        """Format assignment."""
        rhs = self.c_format(expr.rhs)
        lhs = self.c_format(expr.lhs)
        return f"{lhs} {expr.op} {rhs}\n"

    def format_conditional(self, s):
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
        return f"({t} if {c} else {f})"

    def format_symbol(self, s):
        """Format a symbol."""
        return f"{s.name}"

    def format_mathfunction(self, f):
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
        args = [self.c_format(arg) for arg in f.args]
        if "bessel" in function:
            return "0"
        if function == "erf":
            return f"math.erf({args[0]})"
        argstr = ", ".join(args)
        return f"np.{function}({argstr})"

    c_impl = {
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

    def c_format(self, s):
        """Format output."""
        name = s.__class__.__name__
        try:
            return self.c_impl[name](self, s)
        except KeyError:
            raise RuntimeError("Unknown statement: ", name)
