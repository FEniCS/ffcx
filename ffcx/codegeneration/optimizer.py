from typing import List, Union
import ffcx.codegeneration.lnodes as L
from collections import defaultdict
from ffcx.ir.representationutils import QuadratureRule


def optimize(code: List[L.LNode], quadrature_rule) -> List[L.LNode]:
    """Optimize code.

    Parameters
    ----------
    code : list of LNodes
        List of LNodes to optimize.

    Returns
    -------
    list of LNodes
        Optimized list of LNodes.

    """
    # Fuse sections with the same name and same annotations
    code = fuse_sections(code, "Coefficient")
    code = fuse_sections(code, "Jacobian")
    for i, section in enumerate(code):
        if isinstance(section, L.Section):
            if L.Annotation.fuse in section.annotations:
                section = fuse_loops(section)
                code[i] = section
            if L.Annotation.licm in section.annotations:
                section = licm(section, quadrature_rule)
                code[i] = section

    return code


def fuse_sections(code: List[L.LNode], name) -> List[L.LNode]:
    """Fuse sections with the same name and same annotations.

    Parameters
    ----------
    code : list of LNodes
        List of LNodes to fuse.

    Returns
    -------
    list of LNodes
        Fused list of LNodes.

    """
    statements = []
    indices = []
    input = []
    output = []
    annotations = []

    for i, section in enumerate(code):
        if isinstance(section, L.Section):
            if section.name == name:
                statements += section.statements
                indices.append(i)
                input += section.input
                output += section.output
                annotations = section.annotations

    # Remove duplicated inputs
    input = list(set(input))
    # Remove duplicated outputs
    output = list(set(output))

    section = L.Section(name, statements, input, output, annotations)

    # Replace the first section with the fused section
    code = code.copy()
    if indices:
        code[indices[0]] = section
        # Remove the other sections
        code = [c for i, c in enumerate(code) if i not in indices[1:]]

    return code


def fuse_loops(code: L.Section) -> L.Section:
    """Fuse loops with the same range and same annotations.

    Parameters
    ----------
    code : list of LNodes
        List of LNodes to fuse.

    Returns
    -------
    list of LNodes
        Fused list of LNodes.

    """
    loops = defaultdict(list)
    pre_loop = []
    for statement in code.statements:
        if isinstance(statement, L.ForRange):
            id = (statement.index, statement.begin, statement.end)
            loops[id].append(statement.body)
        else:
            pre_loop.append(statement)

    for info, body in loops.items():
        index, begin, end = info
        loop = L.ForRange(index, begin, end, body)
        pre_loop.append(loop)

    return L.Section(code.name, pre_loop, code.input, code.output)


def get_statements(statement: Union[L.Statement, L.StatementList]) -> List[L.LNode]:
    """Get statements from a statement list.

    Parameters
    ----------
    statement : LNode
        Statement list.

    Returns
    -------
    list of LNodes
        List of statements.

    """
    if isinstance(statement, L.StatementList):
        return [statement.expr for statement in statement.statements]
    else:
        return [statement.expr]


def check_dependency(statement: L.Statement, index: L.Symbol) -> bool:
    """Check if a statement depends on a given index.

    Parameters
    ----------
    statement : LNode
        Statement to check.
    index : L.Symbol
        Index to check.

    Returns
    -------
    bool
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
    else:
        raise NotImplementedError(f"Statement {statement} not supported.")

    return False


def licm(section: L.Section, quadrature_rule: QuadratureRule) -> L.Section:
    """Perform loop invariant code motion.

    Parameters
    ----------
    code : list of LNodes
        List of LNodes to optimize.

    Returns
    -------
    list of LNodes
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

    pre_loop: List[L.LNode] = []
    for lhs, rhs in expressions.items():
        for r in rhs:
            hoist_candidates = []
            for arg in r.args:
                dependency = check_dependency(arg, inner_loop.index)
                if not dependency:
                    hoist_candidates.append(arg)
            if (len(hoist_candidates) > 1):
                name = f"temp_{quadrature_rule.id()}_{counter}"
                counter += 1
                temp = L.Symbol(name, L.DataType.SCALAR)
                for h in hoist_candidates:
                    r.args.remove(h)
                # update expression with new temp
                r.args.append(L.ArrayAccess(temp, [outer_loop.index]))
                # create code for hoisted term
                size = outer_loop.end.value - outer_loop.begin.value
                pre_loop.append(L.ArrayDecl(temp, size, [0]))
                body = L.Assign(L.ArrayAccess(temp, [outer_loop.index]), L.Product(hoist_candidates))
                pre_loop.append(L.ForRange(outer_loop.index, outer_loop.begin, outer_loop.end, [body]))

    section.statements = pre_loop + section.statements

    return section
