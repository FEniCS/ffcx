"Code generation for evaluate_dof."

__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2009"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-07

from ffc.codesnippets import jacobian, evaluate_f
from ffc.cpp import format

# Note: Heavily relies on remove_unused to do its job

def _evaluate_dof(ir):
    "Generate code for evaluate_dof"

    # Extract degrees of freedom, function mappings and offsets.
    dofs = ir["dofs"]
    mappings = ir["mappings"]
    offsets = ir["offsets"]

    # Define possibly necessary geometry information
    value_dim = ir["value_dim"]
    cell_dim = ir["cell_dimension"]
    code = jacobian[cell_dim] % {"restriction": ""}

    # Declare variable for physical coordinates
    code += format["comment"]("Declare variable for physical coordinates")
    code += format["declaration"]("double", "y[%d]" % cell_dim)

    # Declare variables for storing the result
    code += format["comment"]("Declare variables for result of evaluation")
    code += format["declaration"]("double", "values[%d]" % value_dim)
    code += format["declaration"]("double", "copy[%d]" % value_dim)
    code += format["declaration"]("double", "result", 0.0)

    # Generate switch bodies for each degree of freedom
    cases = [_generate_body(dof, mappings[i], cell_dim, value_dim, offsets[i])
             for (i, dof) in enumerate(dofs)]

    # Define switch
    code += format["comment"]("Evaluate degree of freedom of function")
    code += format["switch"]("i", cases)

    # Return code as string
    return code

def _generate_body(dof, mapping, cell_dim, value_dim, offset=0):
    "Generate code for a certain dof"

    code = ""

    # Prefetch formats
    iadd = format["iadd"]
    inner_product = format["inner product"]
    component = format["component"]
    assign = format["assign"]
    precision = format["float"]

    # Aid mapping points from reference to physical element
    coefficients = affine_weights(cell_dim)

    # Iterate over the points for this dof. (Integral moments may have several.)
    for point in dof.keys():

        # Map point onto physical element: y = F(x)
        w = coefficients(point)

        for j in range(cell_dim):
            # Compute contraction of weights and coordinates
            value = inner_product(w, [component("x", (k, j))
                                      for k in range(cell_dim + 1)])

            # Assign component j of physical coordinate
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
                values += [component("values" , tokens[1] + offset)]

        # Multiply by weights and sum
        weights = [c[0] for c in dof[point]]
        value = inner_product(weights, values)

        # Add value to result
        code += iadd("result", value)

    # Return result
    code += format["return"]("result")


    return code

def _generate_mapping(mapping, dim, offset):
    "Generate code for mapping function values according to 'mapping'"

    multiply = format["multiply"]
    inner_product = format["inner product"]
    assign = format["assign"]
    component = format["component"]

    code = ""

    if mapping == "affine":
        return code
    elif mapping == "contravariant piola":

        code += format["comment"]("Take inverse of contravariant Piola")

        # Copy of appropriate values
        copies = [assign(component("copy", i), component("values", i+offset))
                  for i in range(dim)]
        code += "".join(copies)

        for i in range(dim):

            # Compute inverse piola for this row
            inv_jacobian_row = ["Jinv_%d%d" % (i, j) for j in range(dim)]
            copies = [component("copy", j) for j in range(dim)]
            values = multiply("detJ", inner_product(inv_jacobian_row, copies))

            # Assign values
            code += assign(component("values", i+offset), values)

    elif mapping == "covariant piola":

        code += format["comment"]("Take inverse of covariant Piola")

        # Copy of appropriate values
        copies = [assign(component("copy", i), component("values", i+offset))
                  for i in range(dim)]
        code += "".join(copies)

        for i in range(dim):
            jacobian_row = ["J_%d%d" % (j, i) for j in range(dim)]
            copies = [component("copy", j) for j in range(dim)]
            values = inner_product(jacobian_row, copies)

            # Assign values
            code += assign(component("values", i+offset), values)

    else:
        raise Exception

    return code

def affine_weights(dim):
    "Compute coefficents for mapping from reference to physical element"

    if dim == 1:
        return lambda x: (1.0 - x[0], x[0])
    elif dim == 2:
        return lambda x: (1.0 - x[0] - x[1], x[0], x[1])
    elif dim == 3:
        return lambda x: (1.0 - x[0] - x[1] - x[2], x[0], x[1], x[2])


def _evaluate_dofs(ir):
    return "// NotImplementedYet"


