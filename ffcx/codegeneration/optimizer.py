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
    """Hashable identity of an expression, used to recognise equal hoisted factors.

    LNodes expressions are not generally hashable, so build a structural key
    that compares equal exactly when two expressions are the same expression.

    Args:
        expr: Expression to key.

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
    # Assignments for every hoisted temp are collected here and emitted as
    # one shared loop below, rather than one separate same-range ForRange
    # per temp: a form with many block entries (e.g. a vector/tensor
    # equation) can hoist dozens of temps, and dozens of tiny loops back to
    # back give the compiler's vectoriser nothing to work with, where one
    # loop computing all of them per iteration is a much better target.
    hoist_loop_body: list[L.LNode] = []
    size = outer_loop.end.value - outer_loop.begin.value
    # Structural key of an already-hoisted expression -> its temp symbol. A
    # vector-valued form with identical, uncoupled per-component blocks
    # (e.g. vector Poisson's diagonal Laplacian blocks) produces the exact
    # same hoisted expression for each component; without this it was
    # recomputed and stored again from scratch for every one of them.
    seen_hoisted: dict[tuple, L.Symbol] = {}
    for lhs, rhs in expressions.items():
        for r in rhs:
            hoist_candidates = []
            for arg in r.args:
                dependency = check_dependency(arg, inner_loop.index)
                if not dependency:
                    hoist_candidates.append(arg)
            if len(hoist_candidates) > 1:
                for h in hoist_candidates:
                    r.args.remove(h)
                hoisted_product = L.Product(hoist_candidates)
                hoisted_key = _key(hoisted_product)
                temp = seen_hoisted.get(hoisted_key)
                if temp is None:
                    # create new temp
                    name = f"temp_{counter}"
                    counter += 1
                    temp = L.Symbol(name, L.DataType.SCALAR)
                    # Every entry gets written unconditionally by the loop
                    # below (a plain Assign, not AssignAdd), so a zero-fill
                    # here is dead.
                    pre_loop.append(L.ArrayDecl(temp, size))
                    hoist_loop_body.append(
                        L.Assign(L.ArrayAccess(temp, [outer_loop.index]), hoisted_product)
                    )
                    seen_hoisted[hoisted_key] = temp
                # update expression with new temp
                r.args.append(L.ArrayAccess(temp, [outer_loop.index]))

    if hoist_loop_body:
        pre_loop.append(
            L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, hoist_loop_body)
        )

    section.statements = pre_loop + section.statements

    return section
