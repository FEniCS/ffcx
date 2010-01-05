"Code generation for evaluate_dof."

__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2009"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-05

from ffc.codesnippets import jacobian
from ffc.cpp import format

# Note: Needs cpp clean-up
# Note: Heavily relies on remove_unused to do its job

def _evaluate_dof(ir):
    "Generate code for evaluate_dof"

    code = []

    # Extract degrees of freedom, function mappings and offsets.
    dofs = ir["dofs"]
    mappings = ir["mappings"]
    offsets = ir["offsets"]

    # Define possibly necessary geometry information
    value_dim = ir["value_dim"]
    cell_dim = ir["cell_dimension"]
    code += [jacobian[cell_dim] % {"restriction": ""}]

    # Declare variable for physical coordinates
    code += ["\n" + format["comment"]("Declare variable for physical coordinates")]
    code += ["double y[%d];" % cell_dim]

    # Declare variables for storing the result
    code += ["\n" + format["comment"]("Declare variables for result of evaluation")]
    code += ["double values[%d];" % value_dim]
    code += ["double copy[%d];" % value_dim]
    code += ["double result = 0;"]

    # Generate switch bodies for each degree of freedom
    cases = [_generate_body(dof, mappings[i], cell_dim, value_dim, offsets[i])
             for (i, dof) in enumerate(dofs)]

    # Define switch
    code += ["\n" + format["comment"]("Evaluate degree of freedom of function")]
    code += [format["switch"]("i", cases)]

    # Return code as string
    return "\n".join(code)

def _generate_body(dof, mapping, cell_dim, value_dim, offset=0):
    "Generate code for a certain dof"

    code = []

    # Prefetch formats
    add = format["add"]
    iadd = format["iadd"]
    multiply = format["multiply"]
    precision = format["float"]

    # Aid mapping points from reference to physical element
    coefficients = affine_weights(cell_dim)

    # Iterate over the points for this dof. (Integral moments may have several.)
    for point in dof.keys():

        # Map point onto physical element: y = F(x)
        w = coefficients(point)
        for j in range(cell_dim):
            value = add([multiply([precision(w[k]), "x[%d][%d]" % (k, j)])
                         for k in range(cell_dim + 1)])
            code += ["y[%d] = %s;" % (j, value)]

        # Evaluate function at physical point
        code += ["f.evaluate(values, y, c);"]

        # Map function values to the reference element
        code += [_generate_mapping(mapping, cell_dim, offset)]

        # Take directional components
        values = []
        for (i, tokens) in enumerate(dof[point]):
            if tokens[1] == ():
                values += ["values[0]"]
            else:
                values += ["values[%d]" % tokens[1]]

        # Multiply by weights and sum
        weights = [precision(c[0]) for c in dof[point]]
        value = add([multiply([weights[i], values[i]]) for i in range(len(values))])

        # Add value to result
        code += [iadd("result", value)]

    # Return result
    code += ["return result;"]

    return "\n".join(code)

def _generate_mapping(mapping, dim, offset):
    "Generate code for mapping function values according to 'mapping'"

    add = format["add"]
    multiply = format["multiply"]

    code = []
    if mapping == "affine":
        pass
    elif mapping == "contravariant piola":
        code += [format["comment"]("Take inverse of contravariant Piola")]
        code += ["copy[%d] = values[%d];" % (i, i+offset) for i in range(dim)]
        values = ["detJ*(" + add([multiply(["Jinv_%d%d" % (i, j), "copy[%d]" % j])
                                  for j in range(dim)]) + ")"
                  for i in range(dim)]
        code += ["values[%d] = %s;" % (i + offset, values[i]) for i in range(dim)]

    elif mapping == "covariant piola":
        code += [format["comment"]("Take inverse of covariant Piola")]
        code += ["copy[%d] = values[%d]" % (i, i+offset) for i in range(dim)]

        values = [add([multiply(["J_%d%d" % (j, i), "copy[%d]" % j])
                       for j in range(dim)]) + ")"
                  for i in range(dim)]
        code += ["values[%d] = %s;" % (i + offset, values[i]) for i in range(dim)]
    else:
        raise Exception

    return "\n".join(code)

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


