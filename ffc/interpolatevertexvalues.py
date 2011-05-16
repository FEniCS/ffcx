"Code generation for interpolate_vertex_values."

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC.  If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Kristian B. Oelgaard 2010
#
# Last changed: 2010-02-09

from ffc.cpp import format, remove_unused

# Extract code manipulation formats
inner =     format["inner product"]
component = format["component"]
assign =    format["assign"]
multiply =  format["multiply"]

# Extract formats for the Jacobians
J =       format["J"]
Jinv =    format["inv(J)"]
invdetJ = format["inverse"](format["det(J)"](None))

f_dof_values =    format["argument dof values"]
f_vertex_values = format["argument vertex values"]

def interpolate_vertex_values(ir):
    "Generate code for interpolate_vertex_values."

    # Handle unsupported elements.
    if isinstance(ir, str):
        return format["exception"]("interpolate_vertex_values: %s" % ir)

    # Add code for Jacobian if necessary
    code = []
    dim = ir["cell_dim"]
    if ir["needs_jacobian"]:
        code.append(format["jacobian and inverse"](dim))

    # Compute total value dimension for (mixed) element
    total_dim = sum(data["value_size"] for data in ir["element_data"])

    # Generate code for each element
    value_offset = 0
    space_offset = 0
    for data in ir["element_data"]:
        # Add vertex interpolation for this element
        code.append(format["comment"]("Evaluate function and change variables"))
        code.append(_interpolate_vertex_values_element(data, dim, total_dim,
                                                       value_offset,
                                                       space_offset))

        # Update offsets for value- and space dimension
        value_offset += data["value_size"]
        space_offset += data["space_dim"]

    # Remove unused variables. (Not tracking set of used variables in
    # order to keep this code clean. Since generated code is of
    # limited size, this should be ok.)
    clean_code = remove_unused("\n".join(code))
    return clean_code


def _interpolate_vertex_values_element(data, dim, total_value_size,
                                       value_offset=0, space_offset=0):

    # Extract vertex values for all basis functions
    vertex_values = data["basis_values"]
    value_size = data["value_size"]
    space_dim = data["space_dim"]
    mapping = data["mapping"]

    # Map basis values according to element mapping. Assumes single
    # mapping for each (non-mixed) element
    change_of_variables = _change_variables(data["mapping"], dim, space_dim)

    # Create code for each value dimension:
    code = []
    for k in range(value_size):

        # Create code for each vertex x_j
        for (j, values_at_vertex) in enumerate(vertex_values):

            if value_size == 1: values_at_vertex = [values_at_vertex]

            # Map basis functions using appropriate mapping
            components = change_of_variables(values_at_vertex, k)

            # Contract coefficients and basis functions
            dof_values = [component(f_dof_values, i + space_offset)
                          for i in range(space_dim)]

            value = inner(dof_values, components)

            # Assign value to correct vertex
            index = j*total_value_size + (k + value_offset)
            code.append(assign(component(f_vertex_values, index), value))

    return "\n".join(code)


def _change_variables(mapping, dim, space_dim):
    """

    How to map a field G from the reference domain to a physical
    domain: For the converse approach -- see evaluatedof.py

    Let g be a field defined on the reference domain domain T_0 with
    reference coordinates X. Let T be a a physical domain with
    coordinates x. Assume that F: T_0 -> T such that

      x = F(X)

    Let J be the Jacobian of F, i.e J = dx/dX and let K denote the
    inverse of the Jacobian K = J^{-1}. Then we (currently) have the
    following three types of mappings:

    'affine' mapping for G:

      g(x) = G(X)

    For vector fields G:

    'contravariant piola' mapping for f:

      g(x) = 1.0/det(J) J G(X)   i.e   g_i(x) = 1.0/det(J) J_ij G_j(X)

    'covariant piola' mapping for f:

      g(x) = K^T G(X)              i.e   g_i(x) = K^T_ij G_j(X) = K_ji G_j(X)
    """

    if mapping is "affine":
        change_of_variables = lambda G, i: G[i]
    elif mapping == "contravariant piola":
        change_of_variables = lambda G, i: [multiply([invdetJ, inner([J(i, j) for j in range(dim)],
                                                                     [G[j][index] for j in range(dim)])])
                                            for index in range(space_dim)]
    elif mapping == "covariant piola":
        change_of_variables = lambda G, i: [inner([Jinv(j, i) for j in range(dim)],
                                                  [G[j][index] for j in range(dim)])
                                            for index in range(space_dim)]
    else:
        raise Exception, "No such mapping: %s accepted" % mapping
    return change_of_variables
