"""Code generation for evaluate_dof.

This module generates the functions evaluate_dof and evaluate_dofs.
These evaluate the degree of freedom (dof) number i and all degrees of
freedom for an element respectively.

Each dof L is assumed to act on a field f in the following manner:

  L(f) = w_{j, k} f_k(x_j)

where w is a set of weights, j is an index set corresponding to the
number of points involved in the evaluation of the functional, and k
is a multi-index set with rank corresponding to the value rank of the
function f.

For common degrees of freedom such as point evaluations and
directional component evaluations, there is just one point. However,
for various integral moments, the integrals are evaluated using
quadrature. The number of points therefore correspond to the
quadrature points.

The points x_j, weights w_{j, k} and components k are extracted from
FIAT (functional.pt_dict) in the intermediate representation stage.

"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2009"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-13

from ffc.codesnippets import jacobian, evaluate_f, cell_coordinates
from ffc.cpp import format

__all__ = ["evaluate_dof", "evaluate_dofs", "affine_weights"]

# Prefetch formats:
comment = format["comment"]
declaration = format["declaration"]
assign = format["assign"]
component = format["component"]
iadd = format["iadd"]
inner = format["inner product"]
add = format["add"]
multiply = format["multiply"]
J = format["J"]
Jinv = format["Jinv"]
detJ = format["det(J)"]

def evaluate_dof(ir):
    "Generate code for evaluate_dof"

    # Define necessary geometry information based on the ir
    code  = _required_declarations(ir)

    # Extract variables
    mappings = ir["mappings"]
    offsets = ir["offsets"]
    cell_dim = ir["cell_dimension"]

    # Generate switch bodies for each degree of freedom
    cases = [_generate_body(dof, mappings[i], cell_dim, offsets[i])
             for (i, dof) in enumerate(ir["dofs"])]

    # Combine with return statements:
    cases = [c + format["return"](r) for (c, r) in cases]

    # Define switch
    code += comment("Evaluate degree of freedom of function")
    code += format["switch"]("i", cases)

    return code


def evaluate_dofs(ir):
    "Generate code for evaluate_dofs"

    # Define necessary geometry information based on the ir
    code = _required_declarations(ir)

    # Extract variables
    mappings = ir["mappings"]
    offsets = ir["offsets"]
    cell_dim = ir["cell_dimension"]

    # Generate code for each degree of freedom
    cases = [_generate_body(dof, mappings[i], cell_dim, offsets[i])
             for (i, dof) in enumerate(ir["dofs"])]

    # Combine with assignments:
    cases = [c + format["assign"]("values[%d]" % i, r)
             for (i, (c, r)) in enumerate(cases)]

    code += "\n".join(cases)
    return code


def _required_declarations(ir):
    """Generate code for declaring required variables and geometry
    information.
    """

    # Declare variable for storing the result
    code = comment("Declare variables for result of evaluation")
    code += declaration("double", "vals[%d]" % ir["value_dim"])

    # Declare variable for physical coordinates
    cell_dim = ir["cell_dimension"]
    code += comment("Declare variable for physical coordinates")
    code += declaration("double", "y[%d]" % cell_dim)

    # Add Jacobian and extra result variable if needed
    needs_jacobian = any(["piola" in m for m in ir["mappings"]])
    if needs_jacobian:
        code += jacobian[cell_dim] % {"restriction": ""}
        code += declaration("double", "result")
    else:
        code += cell_coordinates

    return code


def _generate_body(dof, mapping, cell_dim, offset=0, result="result"):
    "Generate code for a single dof."

    code = ""

    # Make adjustments if multiple points
    points = dof.keys()
    multiple_points = len(points) > 1
    code += assign(result, 0.0) if multiple_points else ""
    compute_result = iadd if multiple_points else assign

    # Aid mapping points from reference to physical element
    coefficients = affine_weights(cell_dim)

    # Iterate over the points for this dof. (Integral moments may have several.)
    for x in points:

        # Map point onto physical element: y = F_K(x)
        w = coefficients(x)
        for j in range(cell_dim):
            y = inner(w, [component("x", (k, j)) for k in range(cell_dim + 1)])
            code += assign(component("y", j), y)

        # Evaluate function at physical point
        code += evaluate_f

        # Map function values to the reference element
        F = _change_variables(mapping, cell_dim, offset)

        # Simple affine functions deserve special case:
        if len(F) == 1 and not multiple_points:
            return (code, F[0])

        # Take inner product between components and weights
        value = add([multiply([w, F[k[0]] if not k == () else F[0]])
                     for (w, k) in dof[x]])

        # Add value to result variable
        code += compute_result(result, value)

    # Return result
    return (code, result)

def _change_variables(mapping, dim, offset):
    """Generate code for mapping function values according to
    'mapping' and offset.
    """

    # meg: Various mappings must be handled both here and in
    # interpolate_vertex_values. Could this be abstracted out?

    if mapping == "affine":
        return [component("vals", offset)]

    elif mapping == "contravariant piola":
        # Map each component from physical to reference using inverse
        # contravariant piola
        values = []
        for i in range(dim):
            inv_jacobian_row = [Jinv(i, j) for j in range(dim)]
            components = [component("vals", j + offset) for j in range(dim)]
            values += [multiply([detJ, inner(inv_jacobian_row, components)])]
        return values

    elif mapping == "covariant piola":
        # Map each component from physical to reference using inverse
        # covariant piola
        values = []
        for i in range(dim):
            jacobian_column = [J(j, i) for j in range(dim)]
            components = [component("vals", j + offset) for j in range(dim)]
            values += [inner(jacobian_column, components)]
        return values
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



