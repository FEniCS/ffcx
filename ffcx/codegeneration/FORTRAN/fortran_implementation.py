
from ffcx.codegeneration.utils import scalar_to_value_type
import ffcx.codegeneration.lnodes as L


class FortranFormatter(object):
    def __init__(self, scalar):
        self.declarations = ""
        self.scalar_type = scalar
        self.real_type = scalar_to_value_type(scalar)

    def format_statement_list(self, slist):
        code = ""
        for s in slist.statements:
            code += self.c_format(s)
        return code

    def format_comment(self, c):
        return "! " + c.comment + "\n"

    def format_array_decl(self, arr):
        types = {"double": "DOUBLE PRECISION", "float": "REAL",
                 "double _Complex": "DOUBLE COMPLEX", "float _Complex": "COMPLEX"}
        t = "VOID"
        if arr.symbol.dtype == L.DataType.SCALAR:
            t = types[self.scalar_type]
        elif arr.symbol.dtype == L.DataType.REAL:
            t = types[self.real_type]

        symbol = self.c_format(arr.symbol)
        dims = ",".join([f"0:{i-1}" for i in arr.sizes])
        decl = f"      {t}, DIMENSION({dims}) :: {symbol}\n"
        self.declarations += decl

        import textwrap
        if arr.values is None:
            return ""
        vals = ", ".join(f"{val}" for val in arr.values.flatten())
        vals = "&\n".join(textwrap.wrap(vals, 75))
        code = f"      {symbol} = reshape((/ &\n{vals}&\n /), shape({symbol}))\n"
        return code

    def format_array_access(self, arr):
        symbol = self.c_format(arr.array)
        idx = ", ".join(self.c_format(ix) for ix in arr.indices)
        return f"{symbol}({idx})"

    def format_multi_index(self, index):
        return self.c_format(index.global_index)

    def format_variable_decl(self, v):
        sym = self.c_format(v.symbol)
        val = self.c_format(v.value)
        return f"      {sym} = {val}\n"

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
        end = self.c_format(r.end - 1)
        index = self.c_format(r.index)
        output = f"      DO {index} = {begin}, {end}\n"
        b = self.c_format(r.body).split("\n")
        output += "  " + "\n  ".join(b)
        output += "    END DO\n"
        return output

    def format_statement(self, s):
        return self.c_format(s.expr)

    def format_assign(self, expr):
        rhs = self.c_format(expr.rhs)
        lhs = self.c_format(expr.lhs)
        if expr.op == "+=":
            return f"      {lhs} = {lhs} + {rhs}\n"
        else:
            return f"      {lhs} {expr.op} {rhs}\n"

    def format_conditional(self, c):
        return "conditional()"

    def format_symbol(self, s):
        return f"{s.name}"

    def format_math_function(self, c):
        lookup = {"fabs": "DABS"}
        fun = lookup.get(c.function, c.function)
        args = ",".join(self.c_format(arg) for arg in c.args)
        return f"{fun}({args})"

    c_impl = {
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
        "MathFunction": format_math_function,
    }

    def c_format(self, s):
        name = s.__class__.__name__
        try:
            st = self.c_impl[name](self, s)
            print(name, st)
            return st
        except KeyError:
            raise RuntimeError("Unknown statement: ", name)
