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

# Copyright (C) 2009 Marie E. Rognes
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Kristian B. Oelgaard 2010-2011
# Modified by Anders Logg 2013
#
# First added:  2009-xx-yy
# Last changed: 2013-01-10

from ffc.cpp import format, remove_unused
from ffc.utils import pick_first

__all__ = ["evaluate_dof_and_dofs", "affine_weights"]

# Prefetch formats:
comment =   format["comment"]
declare =   format["declaration"]
assign =    format["assign"]
component = format["component"]
iadd =      format["iadd"]
inner =     format["inner product"]
add =       format["addition"]
multiply =  format["multiply"]
J =         format["J"]
Jinv =      format["inv(J)"]
detJ =      format["det(J)"](None)
ret =       format["return"]
f_i =       format["argument dof num"]
f_values =  format["argument values"]
f_double =  format["float declaration"]
f_vals =    format["dof vals"]
f_result =  format["dof result"]
f_y =       format["dof physical coordinates"]
f_x =       format["vertex_coordinates"]
f_int =     format["int declaration"]
f_X =       format["dof X"]
f_D =       format["dof D"]
f_W =       format["dof W"]
f_copy =    format["dof copy"]
f_r, f_s =  format["free indices"][:2]
f_loop =    format["generate loop"]

map_onto_physical = format["map onto physical"]

def evaluate_dof_and_dofs(ir):
    "Generate code for evaluate_dof and evaluate_dof."

    # Generate common code
    (reqs, cases) = _generate_common_code(ir)

    # Combine each case with returns for evaluate_dof and switch
    dof_cases = ["%s\n%s" % (c, ret(r)) for (c, r) in cases]
    dof_code = reqs + format["switch"](f_i, dof_cases, ret(format["float"](0.0)))

    # Combine each case with assignments for evaluate_dofs
    dofs_cases = "\n".join("%s\n%s" % (c, format["assign"](component(f_values, i), r))
                           for (i, (c, r)) in enumerate(cases))
    dofs_code = reqs + dofs_cases

    return (dof_code, dofs_code)

def _generate_common_code(ir):

    # Define necessary geometry information based on the ir
    reqs = _required_declarations(ir)

    # Extract variables
    mappings = ir["mappings"]
    offsets  = ir["physical_offsets"]
    gdim = ir["geometric_dimension"]
    tdim = ir["topological_dimension"]

    # Generate bodies for each degree of freedom
    cases = [_generate_body(i, dof, mappings[i], gdim, tdim, offsets[i])
             for (i, dof) in enumerate(ir["dofs"])]

    return (reqs, cases)

def _required_declarations(ir):
    """Generate code for declaring required variables and geometry
    information.
    """
    code = []
    gdim = ir["geometric_dimension"]
    tdim = ir["topological_dimension"]

    # Declare variable for storing the result and physical coordinates
    code.append(comment("Declare variables for result of evaluation"))
    code.append(declare(f_double, component(f_vals, ir["physical_value_size"])))
    code.append("")
    code.append(comment("Declare variable for physical coordinates"))
    code.append(declare(f_double, component(f_y, gdim)))
    code.append("")

    # Check whether Jacobians are necessary.
    needs_inverse_jacobian = any(["contravariant piola" in m
                                  for m in ir["mappings"]])
    needs_jacobian = any(["covariant piola" in m for m in ir["mappings"]])

    # Check if Jacobians are needed
    if not (needs_jacobian or needs_inverse_jacobian):
        return "\n".join(code)

    # Otherwise declare intermediate result variable
    code.append(declare(f_double, f_result))

    # Add sufficient Jacobian information. Note: same criterion for
    # needing inverse Jacobian as for needing oriented Jacobian
    code.append(format["compute_jacobian"](tdim, gdim))
    if needs_inverse_jacobian:
        code.append("")
        code.append(format["compute_jacobian_inverse"](tdim, gdim))
        code.append("")
        code.append(format["orientation"](tdim, gdim))

    return "\n".join(code)

def _generate_body(i, dof, mapping, gdim, tdim, offset=0, result=f_result):
    "Generate code for a single dof."

    points = dof.keys()

    # Generate different code if multiple points. (Otherwise ffc
    # compile time blows up.)
    if len(points) > 1:
        code = _generate_multiple_points_body(i, dof, mapping, gdim, tdim,
                                              offset, result)
        return (code, result)

    # Get weights for mapping reference point to physical
    x = points[0]
    w = affine_weights(tdim)(x)

    # Map point onto physical element: y = F_K(x)
    code = []
    for j in range(gdim):
        y = inner(w, [component(f_x(), (k*gdim + j,)) for k in range(tdim + 1)])
        code.append(assign(component(f_y, j), y))

    # Evaluate function at physical point
    code.append(format["evaluate function"])

    # Map function values to the reference element
    F = _change_variables(mapping, gdim, tdim, offset)

    # Simple affine functions deserve special case:
    if len(F) == 1:
        return ("\n".join(code), multiply([dof[x][0][0], F[0]]))

    # Take inner product between components and weights
    value = add([multiply([w, F[k[0]]]) for (w, k) in dof[x]])

    # Assign value to result variable
    code.append(assign(result, value))
    return ("\n".join(code), result)


def _generate_multiple_points_body(i, dof, mapping, gdim, tdim,
                                   offset=0, result=f_result):

    "Generate c++ for-loop for multiple points (integral bodies)"

    code = [assign(f_result, 0.0)]
    points = dof.keys()
    n = len(points)

    # Get number of tokens per point
    tokens = [dof[x] for x in points]
    len_tokens = pick_first([len(t) for t in tokens])

    # Declare points
    points = format["list"]([format["list"](x) for x in points])
    code += [declare(f_double, component(f_X(i), [n, tdim]),
                     points)]

    # Declare components
    components = [[c[0] for (w, c) in token] for token in tokens]
    components = format["list"]([format["list"](c) for c in components])
    code += [declare(f_int, component(f_D(i), [n, len_tokens]), components)]

    # Declare weights
    weights = [[w for (w, c) in token] for token in tokens]
    weights = format["list"]([format["list"](w) for w in weights])
    code += [declare(f_double, component(f_W(i), [n, len_tokens]), weights)]

    # Declare copy variable:
    code += [declare(f_double, component(f_copy(i), tdim))]

    # Add loop over points
    code += [comment("Loop over points")]

    # Map the points from the reference onto the physical element
    #assert(gdim == tdim), \
    #    "Integral moments not supported for manifolds (yet). Please fix"
    lines_r = [map_onto_physical[tdim][gdim] % {"i": i, "j": f_r}]

    # Evaluate function at physical point
    lines_r.append(comment("Evaluate function at physical point"))
    lines_r.append(format["evaluate function"])

    # Map function values to the reference element
    lines_r.append(comment("Map function to reference element"))
    F = _change_variables(mapping, gdim, tdim, offset)
    lines_r += [assign(component(f_copy(i), k), F_k)
                for (k, F_k) in enumerate(F)]

    # Add loop over directional components
    lines_r.append(comment("Loop over directions"))
    value = multiply([component(f_copy(i),
                                component(f_D(i), (f_r, f_s))),
                      component(f_W(i), (f_r, f_s))])
    # Add value from this point to total result
    lines_s = [iadd(f_result, value)]

    # Generate loop over s and add to r.
    loop_vars_s = [(f_s, 0, len_tokens)]
    lines_r += f_loop(lines_s, loop_vars_s)

    # Generate loop over r and add to code.
    loop_vars_r = [(f_r, 0, n)]
    code += f_loop(lines_r, loop_vars_r)

    code = "\n".join(code)
    return code

def _change_variables(mapping, gdim, tdim, offset):
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
        return [component(f_vals, offset)]

    elif mapping == "contravariant piola":
        # Map each component from physical to reference using inverse
        # contravariant piola
        values = []
        for i in range(tdim):
            inv_jacobian_row = [Jinv(i, j, tdim, gdim) for j in range(gdim)]
            components = [component(f_vals, j + offset) for j in range(gdim)]
            values += [multiply([detJ, inner(inv_jacobian_row, components)])]
        return values

    elif mapping == "covariant piola":
        # Map each component from physical to reference using inverse
        # covariant piola
        values = []
        for i in range(tdim):
            jacobian_column = [J(j, i, gdim, tdim) for j in range(gdim)]
            components = [component(f_vals, j + offset) for j in range(gdim)]
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
