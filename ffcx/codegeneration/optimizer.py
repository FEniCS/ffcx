from typing import List
import ffcx.codegeneration.lnodes as L
from collections import defaultdict


counter = 0


def optimize(code: List[L.LNode]) -> List[L.LNode]:
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
                section = licm(section)
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


def licm(section: L.Section) -> L.Section:
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

    # Check depth of loops
    depth = L.depth(section.statements[0])
    assert depth == 2

    # Get statements in the inner loop
    outer_loop = section.statements[0]
    inner_loop = outer_loop.body.statements[0]
    # get index of inner loops

    index_outer = outer_loop.index
    index_inner = inner_loop.index

    pre_loop: List[L.LNode] = []
    for statement in inner_loop.body.statements:
        if isinstance(statement, L.StatementList):
            continue
        expression = statement.expr
        if isinstance(expression, L.AssignAdd):
            rhs = expression.rhs
            lhs = expression.lhs
            print(hash(lhs))

            # Check if rhs is a sum
            if not isinstance(rhs, L.Sum):
                continue
            for arg in rhs.args:
                hoist = []
                # Check if arg is a product
                if isinstance(arg, L.Product):
                    for factor in arg.args:
                        # Check if factor is ArrayAccess
                        if isinstance(factor, L.ArrayAccess):
                            if index_inner not in factor.indices:
                                hoist.append(factor)
                                # remove from arg.args
                                arg.args.remove(factor)
                        else:
                            hoist.append(factor)
                            arg.args.remove(factor)

                if hoist:
                    # Create new temp
                    temp = L.Symbol(f"t{counter}", L.DataType.REAL)
                    arg.args.append(L.ArrayAccess(temp, [outer_loop.index]))
                    size = outer_loop.end.value - outer_loop.begin.value
                    pre_loop.append(L.ArrayDecl(temp, size, [0]))
                    body = L.Assign(L.ArrayAccess(temp, [outer_loop.index]), L.Product(hoist))
                    pre_loop.append(L.ForRange(index_outer, outer_loop.begin, outer_loop.end, [body]))

    section.statements = pre_loop + section.statements

    return section
