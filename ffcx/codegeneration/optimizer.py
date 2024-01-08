from typing import List
import ffcx.codegeneration.lnodes as L


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

    code = fuse_sections(code, "Jacobian")

    return code


def fuse_sections(code: List[L.LNode], name) -> List[L.LNode]:
    """Fuse sections.

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
    for i, section in enumerate(code):
        if isinstance(section, L.Section):
            if section.name == name:
                statements += section.statements
                indices.append(i)
                input += section.input
                output += section.output

    loops = []
    declarations = []
    for statement in statements:
        if isinstance(statement, L.VariableDecl):
            declarations.append(statement)
        elif isinstance(statement, L.ArrayDecl):
            declarations.append(statement)
        elif isinstance(statement, L.ForRange):
            loops.append(statement)
        else:
            raise NotImplementedError(f"Not expecting {type(statement)}.")

    # Remove duplicated inputs
    input = list(set(input))
    # Remove duplicated outputs
    output = list(set(output))

    section = L.Section(name, declarations + loops, input, output)

    # Replace the first section with the fused section
    code[indices[0]] = section
    # Remove the other sections
    code = [c for i, c in enumerate(code) if i not in indices[1:]]

    return code
