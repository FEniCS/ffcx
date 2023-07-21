import numpy as np


def build_1d_initializer_list(values, formatter, precision=None):
    """Return a list containing a single line formatted like '{ 0.0, 1.0, 2.0 }'."""
    if formatter == str:

        def formatter(x, p):
            return str(x)

    tokens = ["{ "]
    if np.prod(values.shape) > 0:
        sep = ", "
        fvalues = [formatter(v, precision) for v in values]
        for v in fvalues[:-1]:
            tokens.append(v)
            tokens.append(sep)
        tokens.append(fvalues[-1])
    tokens += " }"
    return "".join(tokens)


def build_initializer_lists(values, sizes, level, formatter, precision=None):
    """Return a list of lines with initializer lists for a multidimensional array.

    Example output::

        { { 0.0, 0.1 },
          { 1.0, 1.1 } }

    """
    if formatter == str:

        def formatter(x, p):
            return str(x)

    values = np.asarray(values)
    assert np.prod(values.shape) == np.prod(sizes)
    assert len(sizes) > 0
    assert len(values.shape) > 0
    assert len(sizes) == len(values.shape)
    assert np.all(values.shape == sizes)

    r = len(sizes)
    assert r > 0
    if r == 1:
        return [build_1d_initializer_list(values, formatter, precision=precision)]
    else:
        # Render all sublists
        parts = []
        for val in values:
            sublist = build_initializer_lists(
                val, sizes[1:], level + 1, formatter, precision=precision
            )
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
        for i in range(1, len(lines)):
            lines[i] = "  " + lines[i]
        lines[-1] += " }"
        return lines


def format_statement_list(slist):
    output = ""
    for s in slist.statements:
        output += c_format(s)
    return output


def format_comment(c):
    return "// " + c.comment + "\n"


def format_array_decl(arr):
    symbol = c_format(arr.symbol)
    dims = "".join([f"[{i}]" for i in arr.sizes])
    if arr.values is None:
        vals = "{}"
    else:
        vals = "\n".join(build_initializer_lists(arr.values, arr.sizes, 0, str))
    return f"{arr.typename} {symbol}{dims} = {vals};\n"


def format_array_access(arr):
    name = c_format(arr.array)
    indices = f"[{']['.join(c_format(i) for i in arr.indices)}]"
    return f"{name}{indices}"


def format_variable_decl(v):
    val = c_format(v.value)
    symbol = c_format(v.symbol)
    return f"{v.typename} {symbol} = {val};\n"


def format_nary_op(oper):
    # Format children
    args = [c_format(arg) for arg in oper.args]

    # Apply parentheses
    for i in range(len(args)):
        if oper.args[i].precedence >= oper.precedence:
            args[i] = "(" + args[i] + ")"

    # Return combined string
    return f" {oper.op} ".join(args)


def format_binary_op(oper):
    # Format children
    lhs = c_format(oper.lhs)
    rhs = c_format(oper.rhs)

    # Apply parentheses
    if oper.lhs.precedence >= oper.precedence:
        lhs = f"({lhs})"
    if oper.rhs.precedence >= oper.precedence:
        rhs = f"({rhs})"

    # Return combined string
    return f"{lhs} {oper.op} {rhs}"


def format_literal_float(val):
    return f"{val.value}"


def format_literal_int(val):
    return f"{val.value}"


def format_for_range(r):
    begin = c_format(r.begin)
    end = c_format(r.end)
    index = c_format(r.index)
    output = f"for (int {index} = {begin}; {index} < {end}; ++{index})\n"
    output += "{\n"
    body = c_format(r.body)
    for line in body.split("\n"):
        if len(line) > 0:
            output += f"  {line}\n"
    output += "}\n"
    return output


def format_statement(s):
    return c_format(s.expr)


def format_assign(expr):
    rhs = c_format(expr.rhs)
    lhs = c_format(expr.lhs)
    return f"{lhs} {expr.op} {rhs};\n"


def format_conditional(s):
    # Format children
    c = c_format(s.condition)
    t = c_format(s.true)
    f = c_format(s.false)

    # Apply parentheses
    if s.condition.precedence >= s.precedence:
        c = "(" + c + ")"
    if s.true.precedence >= s.precedence:
        t = "(" + t + ")"
    if s.false.precedence >= s.precedence:
        f = "(" + f + ")"

    # Return combined string
    return c + " ? " + t + " : " + f


def format_symbol(s):
    return f"{s.name}"


def format_math_function(c):
    args = ",".join(c_format(arg) for arg in c.args)
    return f"{c.function}({args})"


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


def c_format(s):
    name = s.__class__.__name__
    try:
        return c_impl[name](s)
    except KeyError:
        raise RuntimeError("Unknown statement: ", name)
