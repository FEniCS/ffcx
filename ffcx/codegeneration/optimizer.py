"""Optimizer."""

from collections import defaultdict

import ffcx.codegeneration.lnodes as L
from ffcx.ir.representationutils import QuadratureRule


def optimize(code: list[L.LNode], quadrature_rule: QuadratureRule) -> list[L.LNode]:
    """Optimize code.

    Args:
        code: List of LNodes to optimize.
        quadrature_rule: TODO.

    Returns:
        Optimized list of LNodes.
    """
    # Fuse sections with the same name and same annotations
    code = fuse_sections(code, "Coefficient")
    code = fuse_sections(code, "Jacobian")
    for i, section in enumerate(code):
        if isinstance(section, L.Section):
            if L.Annotation.fuse in section.annotations:
                section = fuse_loops(section)
            if L.Annotation.licm in section.annotations:
                section = licm(section, quadrature_rule)
            code[i] = section

    return code


def fuse_sections(code: list[L.LNode], name: str) -> list[L.LNode]:
    """Fuse sections with the same name.

    Args:
        code: List of LNodes to fuse.
        name: Common name used by the sections that should be fused

    Returns:
        Fused list of LNodes.
    """
    statements: list[L.LNode] = []
    indices: list[int] = []
    input: list[L.Symbol] = []
    output: list[L.Symbol] = []
    declarations: list[L.Declaration] = []
    annotations: list[L.Annotation] = []

    for i, section in enumerate(code):
        if isinstance(section, L.Section):
            if section.name == name:
                declarations.extend(section.declarations)
                statements.extend(section.statements)
                indices.append(i)
                input.extend(section.input)
                output.extend(section.output)
                annotations = section.annotations

    # Remove duplicated inputs
    input = list(set(input))
    # Remove duplicated outputs
    output = list(set(output))

    section = L.Section(name, statements, declarations, input, output, annotations)

    # Replace the first section with the fused section
    code = code.copy()
    if indices:
        code[indices[0]] = section
        # Remove the other sections
        code = [c for i, c in enumerate(code) if i not in indices[1:]]

    return code


def fuse_loops(code: L.Section) -> L.Section:
    """Fuse loops with the same range and same annotations.

    Args:
        code: List of LNodes to fuse.

    Returns:
        Fused list of LNodes.
    """
    loops = defaultdict(list)
    output_code = []
    for statement in code.statements:
        if isinstance(statement, L.ForRange):
            id = (statement.index, statement.begin, statement.end)
            loops[id].append(statement.body)
        else:
            output_code.append(statement)

    for range, body in loops.items():
        output_code.append(L.ForRange(*range, body))

    return L.Section(code.name, output_code, code.declarations, code.input, code.output)


def get_statements(statement: L.Statement | L.StatementList) -> list[L.LNode]:
    """Get statements from a statement list.

    Args:
        statement: Statement list.

    Returns:
        List of statements.
    """
    if isinstance(statement, L.StatementList):
        return [statement.expr for statement in statement.statements]
    else:
        return [statement.expr]


def check_dependency(statement: L.Statement, index: L.Symbol) -> bool:
    """Check if a statement depends on a given index.

    Args:
        statement: Statement to check.
        index: Index to check.

    Returns:
        True if statement depends on index, False otherwise.
    """
    if isinstance(statement, L.ArrayAccess):
        if index in statement.indices:
            return True
        else:
            for i in statement.indices:
                if isinstance(i, L.Sum) or isinstance(i, L.Product):
                    if index in i.args:
                        return True
    elif isinstance(statement, L.Symbol):
        return False
    elif isinstance(statement, L.LiteralFloat) or isinstance(statement, L.LiteralInt):
        return False
    else:
        raise NotImplementedError(f"Statement {statement} not supported.")

    return False


def _key(expr) -> tuple:
    """Hashable identity of an expression, used to group equal factors.

    LNodes expressions are not generally hashable, so build a structural key
    that compares equal exactly when two factors are the same expression.

    Args:
        expr: Expression appearing as a factor in a product.

    Returns:
        Hashable structural key.
    """
    if isinstance(expr, L.ArrayAccess):
        return ("array", expr.array.name, *(_key(i) for i in expr.indices))
    elif isinstance(expr, L.Symbol):
        return ("symbol", expr.name)
    elif isinstance(expr, L.LiteralFloat | L.LiteralInt):
        return ("literal", expr.value)
    args = getattr(expr, "args", None)
    if args is not None:
        return (type(expr).__name__, *(_key(a) for a in args))
    return (type(expr).__name__, repr(expr))


def licm(section: L.Section, quadrature_rule: QuadratureRule) -> L.Section:
    """Perform loop invariant code motion.

    Args:
        section: List of LNodes to optimize.
        quadrature_rule: TODO.

    Returns:
        Optimized list of LNodes.
    """
    assert L.Annotation.licm in section.annotations

    counter = 0

    # Check depth of loops
    depth = L.depth(section.statements[0])
    if depth != 2:
        return section

    # Get statements in the inner loop
    outer_loop = section.statements[0]
    inner_loop = outer_loop.body.statements[0]

    # Collect all expressions in the inner loop by corresponding RHS
    expressions = defaultdict(list)
    for body in inner_loop.body.statements:
        statements = get_statements(body)
        assert isinstance(statements, list)
        for statement in statements:
            assert isinstance(statement, L.AssignAdd)  # Expecting AssignAdd
            rhs = statement.rhs
            assert isinstance(rhs, L.Product)  # Expecting Sum
            lhs = statement.lhs
            assert isinstance(lhs, L.ArrayAccess)  # Expecting ArrayAccess
            expressions[lhs].append(rhs)

    pre_loop: list[L.LNode] = []
    inner_body: list[L.LNode] = []
    size = outer_loop.end.value - outer_loop.begin.value
    for lhs, rhs in expressions.items():
        # Group the terms of this entry by the factors that stay in the inner
        # loop. Terms sharing them differ only in what is hoisted, so one
        # temporary can hold their sum and the entry needs a single update.
        groups: dict[tuple, tuple[list, list]] = {}
        for r in rhs:
            hoist_candidates = []
            keep = []
            for arg in r.args:
                if check_dependency(arg, inner_loop.index):
                    keep.append(arg)
                else:
                    hoist_candidates.append(arg)
            if len(hoist_candidates) > 1:
                key = tuple(_key(k) for k in keep)
                groups.setdefault(key, (keep, []))[1].append(L.Product(hoist_candidates))
            else:
                # Nothing worth hoisting, emit the term unchanged
                inner_body.append(L.AssignAdd(lhs, r))

        for keep, hoisted in groups.values():
            # Create new temp holding the sum of the hoisted terms
            name = f"temp_{counter}"
            counter += 1
            temp = L.Symbol(name, L.DataType.SCALAR)
            pre_loop.append(L.ArrayDecl(temp, size, [0]))
            rhs_hoisted = hoisted[0] if len(hoisted) == 1 else L.Sum(hoisted)
            body = L.Assign(L.ArrayAccess(temp, [outer_loop.index]), rhs_hoisted)
            pre_loop.append(L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, [body]))
            access = L.ArrayAccess(temp, [outer_loop.index])
            inner_body.append(L.AssignAdd(lhs, L.Product([*keep, access])))

    # Rebuild the loop nest around the regrouped body
    new_inner = L.ForRange(inner_loop.index, inner_loop.begin, inner_loop.end, inner_body)
    new_outer = L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, [new_inner])
    section.statements = pre_loop + [new_outer] + section.statements[1:]

    return section
