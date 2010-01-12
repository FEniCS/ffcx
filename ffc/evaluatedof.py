"Code generation for evaluate_dof."

__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2009"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-12

from ffc.codesnippets import jacobian, evaluate_f, cell_coordinates
from ffc.cpp import format

# Prefetch formats
iadd = format["iadd"]
multiply = format["multiply"]
inner_product = format["inner product"]
component = format["component"]
assign = format["assign"]
precision = format["float"]
declaration = format["declaration"]
comment = format["comment"]
J = format["J"]
Jinv = format["Jinv"]
detJ = format["det(J)"]

def _evaluate_dof(ir):
    "Generate code for evaluate_dof"

    # Define necessary geometry information
    code  = _required_declarations(ir)

    # Extract variables
    mappings = ir["mappings"]
    offsets = ir["offsets"]
    cell_dim = ir["cell_dimension"]
    value_dim = ir["value_dim"]

    # Generate switch bodies for each degree of freedom
    cases = [_generate_body(dof, mappings[i], cell_dim, value_dim, offsets[i])
             for (i, dof) in enumerate(ir["dofs"])]

    # Define switch
    code += comment("Evaluate degree of freedom of function")
    code += format["switch"]("i", cases)

    return code


def _evaluate_dofs(ir):
    "Generate code for evaluate_dofs"
    return "// NotImplementedYet"


def _required_declarations(ir):
    """Generate code for declaring required variables and geometry
    information.
    """

    code = ""

    # Declare variables for storing the result
    value_dim = ir["value_dim"]
    code += comment("Declare variables for result of evaluation")
    code += declaration("double", "values[%d]" % value_dim)
    code += declaration("double", "result", 0.0)

    # Declare variable for physical coordinates
    cell_dim = ir["cell_dimension"]
    code += comment("Declare variable for physical coordinates")
    code += declaration("double", "y[%d]" % cell_dim)

    # Add Jacobian if needed
    needs_jacobian = any(["piola" in m for m in ir["mappings"]])
    if needs_jacobian:
        code += jacobian[cell_dim] % {"restriction": ""}
        code += declaration("double", "copy[%d]" % value_dim)
    else:
        code += cell_coordinates
    return code

def _generate_body(dof, mapping, cell_dim, value_dim, offset=0):
    "Generate code for a single dof."

    code = ""

    # Aid mapping points from reference to physical element
    coefficients = affine_weights(cell_dim)

    # Iterate over the points for this dof. (Integral moments may have several.)
    for point in dof.keys():

        # Map point onto physical element: y = F(x)
        w = coefficients(point)
        for j in range(cell_dim):
            value = inner_product(w, [component("x", (k, j))
                                      for k in range(cell_dim + 1)])
            code += assign(component("y", j), value)

        # Evaluate function at physical point
        code += evaluate_f

        # Map function values to the reference element
        code += _generate_mapping(mapping, cell_dim, offset)

        # Take directional components
        values = []
        for (i, tokens) in enumerate(dof[point]):
            if tokens[1] == ():
                values += [component("values", offset)]
            else:
                values += [component("values" , tokens[1][0] + offset)]

        # Take inner product with weights
        weights = [c[0] for c in dof[point]]
        value = inner_product(weights, values)

        # Add value to result variable
        code += iadd("result", value)

    # Return result
    code += format["return"]("result")
    return code

def _generate_mapping(mapping, dim, offset):
    "Generate code for mapping function values according to 'mapping'"

    # meg: Various mappings must be handled both here and in
    # interpolate_vertex_values. Could this be abstracted out?
    code = ""

    if mapping == "affine":
        # Do nothing
        return code

    elif mapping == "contravariant piola":
        # Copy values to be reassigned
        code += "".join([assign(component("copy", i), component("values", i+offset))
                         for i in range(dim)])

        # Map each component from physical to reference using inverse
        # contravariant piola
        code += comment("Take inverse of contravariant Piola")
        for i in range(dim):
            inv_jacobian_row = [Jinv(i, j) for j in range(dim)]
            copies = [component("copy", j) for j in range(dim)]
            values = multiply([detJ, inner_product(inv_jacobian_row, copies)])
            code += assign(component("values", i+offset), values)

    elif mapping == "covariant piola":
        # Copy values to be reassigned
        code += "".join([assign(component("copy", i), component("values", i+offset))
                         for i in range(dim)])

        # Map each component from physical to reference using inverse
        # covariant piola
        code += comment("Take inverse of covariant Piola")
        for i in range(dim):
            jacobian_column = [J(j, i) for j in range(dim)]
            copies = [component("copy", j) for j in range(dim)]
            values = inner_product(jacobian_column, copies)
            code += assign(component("values", i+offset), values)

    else:
        raise Exception, "The mapping (%s) is not allowed" % mapping

    return code

def affine_weights(dim):
    "Compute coefficents for mapping from reference to physical element"

    if dim == 1:
        return lambda x: (1.0 - x[0], x[0])
    elif dim == 2:
        return lambda x: (1.0 - x[0] - x[1], x[0], x[1])
    elif dim == 3:
        return lambda x: (1.0 - x[0] - x[1] - x[2], x[0], x[1], x[2])



