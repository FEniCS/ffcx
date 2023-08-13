import numpy as np

class NumbaFormatter(object):
    def __init__(self, scalar) -> None:
        self.scalar_type = scalar

    def format_statement_list(self, slist):
        output = ""
        for s in slist.statements:
            output += self.c_format(s)
        return output


    def format_comment(self, c):
        return "# " + c.comment + "\n"


    def format_array_decl(self, arr):
        typename = arr.typename
        if "double" in typename:
            dtype = "np.float64"
        elif "float" in typename:
            dtype = "np.float32"
        symbol = self.c_format(arr.symbol)
        if arr.values is None:
            return f"{symbol} = np.empty({arr.sizes}, dtype={dtype})\n"
        av = repr(arr.values)
        if av.startswith("array"):
            av = "np." + av[:-1] + f", dtype={dtype})"
            return f"{symbol} = {av}\n"
        else:
            raise RuntimeError("Not sure about array:", av)


    def format_array_access(self, arr):
        array = self.c_format(arr.array)
        idx = ", ".join(self.c_format(ix) for ix in arr.indices)
        return f"{array}[{idx}]"


    def format_variable_decl(self, v):
        sym = self.c_format(v.symbol)
        val = self.c_format(v.value)
        return f"{sym} = {val}\n"


    def format_nary_op(self, oper):
        # Format children
        args = [self.c_format(arg) for arg in oper.args]

        # Apply parentheses
        for i in range(len(args)):
            if oper.args[i].precedence >= oper.precedence:
                args[i] = "(" + args[i] + ")"

        # Return combined string
        return f" {oper.op} ".join(args)


    def format_binary_op(self, oper):
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
        arg = self.c_format(val.arg)
        return f"-{arg}"


    def format_literal_float(self, val):
        return f"{val.value}"


    def format_literal_int(self, val):
        return f"{val.value}"


    def format_for_range(self, r):
        begin = self.c_format(r.begin)
        end = self.c_format(r.end)
        index = self.c_format(r.index)
        output = f"for {index} in range({begin}, {end}):\n"
        b = self.c_format(r.body).split("\n")
        for line in b:
            output += f"    {line}\n"
        return output


    def format_statement(self, s):
        return self.c_format(s.expr)


    def format_assign(self, expr):
        rhs = self.c_format(expr.rhs)
        lhs = self.c_format(expr.lhs)
        return f"{lhs} {expr.op} {rhs}\n"


    def format_conditional(self, c):
        return "conditional()"


    def format_symbol(self, s):
        return f"{s.name}"


    def format_mathfunction(self, f):
        function = f.function
        args = ", ".join(self.c_format(arg) for arg in f.args)
        return f"{function}({args})"


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
        "Sub": format_binary_op,
        "Mul": format_binary_op,
        "Div": format_binary_op,
        "Neg": format_neg,
        "LiteralFloat": format_literal_float,
        "LiteralInt": format_literal_int,
        "Symbol": format_symbol,
        "Conditional": format_conditional,
        "MathFunction": format_mathfunction,
    }


    def c_format(self, s):
        name = s.__class__.__name__
        try:
            return self.c_impl[name](self, s)
        except KeyError:
            raise RuntimeError("Unknown statement: ", name)
