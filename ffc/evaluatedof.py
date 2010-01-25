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

from ffc.cpp import format, remove_unused, IndentControl
from ffc.utils import pick_first

__all__ = ["evaluate_dof_and_dofs", "affine_weights"]

# Prefetch formats:
comment = format["comment"]
declare = format["declaration"]
assign = format["assign"]
component = format["component"]
iadd = format["iadd"]
inner = format["inner product"]
add = format["addition"]
multiply = format["multiply"]
J = format["J"]
Jinv = format["inv(J)"]
detJ = format["det(J)"]("")
map_onto_physical = format["map onto physical"]

def evaluate_dof_and_dofs(ir):
    "Generate code for evaluate_dof and evaluate_dof."

    # Generate common code
    (reqs, cases) = _generate_common_code(ir)

    # Combine each case with returns for evaluate_dof and switch
    dof_cases = ["%s\n%s" % (c, format["return"](r)) for (c, r) in cases]
    dof_code = reqs + format["switch"]("i", dof_cases, format["return"]("0.0"))

    # Combine each case with assignments for evaluate_dofs
    dofs_cases = "\n".join("%s\n%s" % (c, format["assign"]("values[%d]" % i, r))
                           for (i, (c, r)) in enumerate(cases))
    dofs_code = reqs + dofs_cases

    return (dof_code, dofs_code)

def _generate_common_code(ir):

    # Define necessary geometry information based on the ir
    reqs = _required_declarations(ir)

    # Extract variables
    mappings = ir["mappings"]
    offsets = ir["offsets"]
    cell_dim = ir["cell_dimension"]

    # Generate bodies for each degree of freedom
    cases = [_generate_body(i, dof, mappings[i], cell_dim, offsets[i])
             for (i, dof) in enumerate(ir["dofs"])]

    return (reqs, cases)

def _required_declarations(ir):
    """Generate code for declaring required variables and geometry
    information.
    """
    code = []

    # Declare variable for storing the result and physical coordinates
    cell_dim = ir["cell_dimension"]
    code.append(comment("Declare variables for result of evaluation"))
    code.append(declare("double", "vals[%d]" % ir["value_dim"]))
    code.append(comment("Declare variable for physical coordinates"))
    code.append(declare("double", "y[%d]" % cell_dim))

    # Check whether Jacobians are necessary.
    needs_inverse_jacobian = any(["contravariant piola" in m
                                  for m in ir["mappings"]])
    needs_jacobian = any(["covariant piola" in m for m in ir["mappings"]])

    # Add cell coordinates only if sufficient
    if not (needs_jacobian or needs_inverse_jacobian):
        code.append(format["cell coordinates"])
        return "\n".join(code)

    # Otherwise declare intermediate result variable
    code.append(declare("double", "result"))

    # Add sufficient Jacobian information
    if needs_inverse_jacobian:
        code.append(format["jacobian and inverse"](cell_dim))
    else:
        code.append(format["jacobian"](cell_dim))

    return "\n".join(code)


def _generate_body(i, dof, mapping, cell_dim, offset=0, result="result"):
    "Generate code for a single dof."

    points = dof.keys()

    # Generate different code if multiple points. (Otherwise ffc
    # compile time blows up.)
    if len(points) > 1:
        code = _generate_multiple_points_body(i, dof, mapping,
                                              cell_dim, offset, result)
        return (code, result)

    # Get weights for mapping reference point to physical
    x = points[0]
    w = affine_weights(cell_dim)(x)

    # Map point onto physical element: y = F_K(x)
    code = []
    for j in range(cell_dim):
        y = inner(w, [component("x", (k, j)) for k in range(cell_dim + 1)])
        code.append(assign(component("y", j), y))

    # Evaluate function at physical point
    code.append(format["evaluate function"])

    # Map function values to the reference element
    F = _change_variables(mapping, cell_dim, offset)

    # Simple affine functions deserve special case:
    if len(F) == 1:
        return ("\n".join(code), multiply([dof[x][0][0], F[0]]))

    # Take inner product between components and weights
    value = add([multiply([w, F[k[0]]]) for (w, k) in dof[x]])

    # Assign value to result variable
    code.append(assign(result, value))
    return ("\n".join(code), result)


def _generate_multiple_points_body(i, dof, mapping,
                                   cell_dim, offset=0, result="result"):

    "Generate c++ for-loop for multiple points (integral bodies)"

    code = [assign("result", 0.0)]
    points = dof.keys()
    n = len(points)

    # Get number of tokens per point
    tokens = [dof[x] for x in points]
    len_tokens = pick_first([len(t) for t in tokens])


    # Declare points
    points = format["list"]([format["list"](x) for x in points])
    code += [declare("double", "X_%d[%d][%d]" % (i, n, cell_dim), points)]

    # Declare components
    components = [[c[0] for (w, c) in token] for token in tokens]
    components = format["list"]([format["list"](c) for c in components])
    code += [declare("int", "D_%d[%d][%d]" % (i, n, len_tokens), components)]

    # Declare weights
    weights = [[w for (w, c) in token] for token in tokens]
    weights = format["list"]([format["list"](w) for w in weights])
    code += [declare("double", "W_%d[%d][%d]" % (i, n, len_tokens), weights)]

    # Declare copy variable:
    code += [declare("double", "copy_%d[%d]" % (i, cell_dim))]

    # Add loop over points
    code += [comment("// Loop over points")]
    code += ["for(unsigned int j = 0; j < %d; j++) {\n" % n]

    Indent = IndentControl()
    Indent.increase()

    # Map the points from the reference onto the physical element
    code += [Indent.indent(map_onto_physical[cell_dim] % {"i": i, "j": "j"})]

    # Evaluate function at physical point
    code += [Indent.indent(comment("// Evaluate function at physical point"))]
    code.append(Indent.indent(format["evaluate function"]))

    # Map function values to the reference element
    code += [Indent.indent(comment("// Map function to reference element"))]
    F = _change_variables(mapping, cell_dim, offset)
    code += [Indent.indent(assign(component("copy_%d" % i, k), F_k))
             for (k, F_k) in enumerate(F)]

    # Add loop over directional components
    code += [Indent.indent(comment("// Loop over directions"))]
    code += [Indent.indent("for(unsigned int k = 0; k < %d; k++) {" % len_tokens)]
    Indent.increase()
    value = multiply([component("copy", component("D_%d" % i, ("j", "k"))),
                      component("W_%d" %i, ("j", "k"))])

    # Add value from this point to total result
    code += [Indent.indent(iadd("result", value))]
    code += [Indent.indent("}")]
    Indent.decrease()

    code += [Indent.indent("}")]
    code = "\n".join(code)
    return code


def _change_variables(mapping, dim, offset):
    """Generate code for mapping function values according to
    'mapping' and offset.

    The basics of how to map a field from a physical to the reference
    domain. (For the inverse approach -- see interpolatevertexvalues)

    Let g be a field defined on a physical domain T with physical
    coordinates x. Let T_0 be a reference domain with coordinates
    X. Assume that F: T_0 -> T such that

      x = F(X)

    Let J be the Jacobian of F, i.e J = dx/dX and let K denote the
    inverse of the Jacobian K = J^{-1}. Then we (currently) have the
    following three types of mappings:

    'affine' mapping for g:

      G(X) = g(x)

    For vector fields g:

    'contravariant piola' mapping for g:

      G(X) = det(J) K g(x)   i.e  G_i(X) = det(J) K_ij g_j(x)

    'covariant piola' mapping for g:

      G(X) = J^T g(x)          i.e  G_i(X) = J^T_ij g(x) = J_ji g_j(x)
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



