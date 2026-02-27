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
                annotations = list(section.annotations)

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

    return L.Section(code.name, output_code, code.declarations, list(code.input), list(code.output))


def get_statements(statement: L.Statement | L.StatementList) -> list[L.LExpr]:
    """Get statements from a statement list.

    Args:
        statement: Statement list.

    Returns:
        List of expression nodes.
    """
    if isinstance(statement, L.StatementList):
        return [st.expr for st in statement.statements if isinstance(st, L.Statement)]
    else:
        return [statement.expr]


def check_dependency(statement: L.LExpr, index: L.Symbol | L.MultiIndex) -> bool:
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
    assert isinstance(outer_loop, L.ForRange)
    inner_loop = outer_loop.body.statements[0]
    assert isinstance(inner_loop, L.ForRange)

    # Collect all expressions in the inner loop by corresponding RHS
    expressions: defaultdict[L.ArrayAccess, list[L.Product]] = defaultdict(list)
    for body in inner_loop.body.statements:
        stmts = get_statements(body)
        assert isinstance(stmts, list)
        for stmt in stmts:
            assert isinstance(stmt, L.AssignAdd)  # Expecting AssignAdd
            rhs = stmt.rhs
            assert isinstance(rhs, L.Product)  # Expecting Product
            lhs = stmt.lhs
            assert isinstance(lhs, L.ArrayAccess)  # Expecting ArrayAccess
            expressions[lhs].append(rhs)

    # Build replacement mapping: old Product id → new Product
    # and collect pre-loop declarations/loops
    replacements: dict[int, L.Product] = {}  # id(old_product) → new_product
    pre_loop: list[L.LNode] = []
    for expr_lhs, expr_rhs_list in expressions.items():
        for r in expr_rhs_list:
            hoist_candidates: list[L.LExpr] = []
            for arg in r.args:
                dependency = check_dependency(arg, inner_loop.index)
                if not dependency:
                    hoist_candidates.append(arg)
            if len(hoist_candidates) > 1:
                # create new temp
                name = f"temp_{counter}"
                counter += 1
                temp = L.Symbol(name, L.DataType.SCALAR)
                # Build new Product: keep non-hoisted args + new temp access
                remaining_args: list[L.LExpr] = [a for a in r.args if a not in hoist_candidates]
                remaining_args.append(L.ArrayAccess(temp, [outer_loop.index]))
                replacements[id(r)] = L.Product(tuple(remaining_args))
                # create code for hoisted term
                assert isinstance(outer_loop.end, L.LiteralInt)
                assert isinstance(outer_loop.begin, L.LiteralInt)
                size = outer_loop.end.value - outer_loop.begin.value
                pre_loop.append(L.ArrayDecl(temp, size, [0]))
                hoist_body = L.Assign(
                    L.ArrayAccess(temp, [outer_loop.index]),
                    L.Product(tuple(hoist_candidates)),
                )
                pre_loop.append(
                    L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, [hoist_body])
                )

    if not replacements:
        return section

    # Rebuild inner loop body with replaced Products
    new_inner_stmts: list[L.LNode] = []
    for body in inner_loop.body.statements:
        stmts = get_statements(body)
        new_stmts: list[L.AssignAdd] = []
        for stmt in stmts:
            assert isinstance(stmt, L.AssignAdd)
            new_rhs = replacements.get(id(stmt.rhs), stmt.rhs)
            new_stmt = L.AssignAdd(stmt.lhs, new_rhs)
            new_stmts.append(new_stmt)
        if len(new_stmts) == 1:
            new_inner_stmts.append(new_stmts[0])
        else:
            new_inner_stmts.append(L.StatementList(new_stmts))

    # Rebuild loops from inside out
    new_inner_loop = L.ForRange(inner_loop.index, inner_loop.begin, inner_loop.end, new_inner_stmts)
    new_outer_loop = L.ForRange(
        outer_loop.index, outer_loop.begin, outer_loop.end, [new_inner_loop]
    )

    # Build new section with pre-loop code prepended
    new_statements = pre_loop + [new_outer_loop] + list(section.statements[1:])
    return L.Section(
        section.name,
        new_statements,
        section.declarations,
        section.input,
        section.output,
        section.annotations,
    )
