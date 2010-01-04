""". """

__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2009 "
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-04

from cpp import format

def _evaluate_dof(ir):

    code = []

    # FIXME:
    value_dimension = 1
    cell_dimension = 2

    # Get coordinates defining cell
    code += ["const double * const * x = c.coordinates;"]

    # Declare various variable for physical points
    code += ["double y[%d];" % cell_dimension]
    code += ["double result;"]
    code += ["double values[%d];" % value_dimension]

    # Generate switch bodies for each dof
    switch_bodies = [_generate_body(dof) for dof in ir]

    # Create switch
    code += [format["switch"]("i", switch_bodies)]

    code = "\n".join(code)

    print "code :\n ", code
    return code

def _generate_body(dof):

    code = []

    # Prefetch formats
    add = format["add"]
    multiply = format["multiply"]

    # FIXME
    cell_dimension = 2

    # Iterate over the points (Assume only one for now
    point_index = 0
    point = dof.keys()[point_index]

    # Map reference point to physical point
    w = affine_weights(point)
    for j in range(cell_dimension):
        name = "y[%d]" % j
        value = add([multiply(["%0.9g" % w[k], "x[%d][%d]" % (k, j)])
                     for k in range(cell_dimension + 1)])
        code += [name + " = " + value + ";"]

    # Evaluate function at appropriate point
    code += ["f.evaluate(values, y, c);"]

    # Map function to the reference element using appropriate
    # mapping

    # Take directional components (if appropriate)
    weights, components = dof[point][point_index]
    if components == ():
        value = "values"
    else:
        value = "values[%d]" % component[0]

    value = multiply([value, "%0.9g" % weights])

    # Return result
    if value == "values":
        code += ["return values;"]
    else:
        code += ["result = " + value + ";"]
        code += ["return result;"]

    return "\n".join(code)

def _evaluate_dofs(ir):

    return ""


def affine_weights(x):
    # FIXME: 2D

    return (1.0 - x[0] - x[1], x[0], x[1])

