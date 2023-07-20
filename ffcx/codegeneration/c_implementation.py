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
        output += cs_format(s)
    return output


def format_comment(c):
    return "// " + c.comment + "\n"


def format_array_decl(arr):
    dims = "".join([f"[{i}]" for i in arr.sizes])
    if arr.values is None:
        vals = "{}"
    else:
        vals = "\n".join(build_initializer_lists(arr.values, arr.sizes, 0, str))
    return f"{arr.typename} {arr.symbol}{dims} = {vals};\n"


def format_variable_decl(v):
    return f"{v.typename} {v.symbol} = {v.value};\n"


def format_for_range(r):
    output = f"for (int {r.index} = {r.begin}; {r.index} < {r.end}; ++{r.index})\n"
    output += "{\n"
    output += cs_format(r.body)
    output += "}\n"
    return output


def format_statement(s):
    return cs_format(s.expr)


def format_assign(expr):
    return f"{expr.lhs} {expr.op} {expr.rhs};\n"


c_impl = {
    "StatementList": format_statement_list,
    "Comment": format_comment,
    "ArrayDecl": format_array_decl,
    "VariableDecl": format_variable_decl,
    "ForRange": format_for_range,
    "Statement": format_statement,
    "Assign": format_assign,
    "AssignAdd": format_assign,
}


def cs_format(s):
    name = s.__class__.__name__
    try:
        return c_impl[name](s)
    except KeyError:
        raise RuntimeError("Unknown statement: ", type(s))
