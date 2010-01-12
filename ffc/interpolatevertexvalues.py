"Code generation for interpolate_vertex_values."

__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2009"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-12

from ffc.codesnippets import jacobian
from ffc.cpp import format

# Extract code manipulation formats
inner_product = format["inner product"]
component = format["component"]
precision = format["float"]
assign = format["assign"]
multiply = format["multiply"]

# Extract formats for the Jacobians
J = format["J"]
Jinv = format["Jinv"]
invdetJ = "1.0/%s" % format["det(J)"]

def _interpolate_vertex_values(ir):
    "Generate code for interpolate_vertex_values."

    # meg: Note: I think this works, but it might be more elegant to
    # give mixed element a tabulate function

    code = ""

    # Add code for Jacobian if necessary
    dim = ir["cell_dim"]
    if ir["needs_jacobian"]:
        code += jacobian[dim] % {"restriction": ""}

    # Compute total value dimension for (mixed) element
    total_value_dim = sum(data["value_dim"] for data in ir["element_data"])

    # Generate code for each element
    value_offset = 0
    space_offset = 0
    for data in ir["element_data"]:

        code += format["comment"]("Evaluate at vertices and map")
        code += _interpolate_vertex_values_element(data, dim, total_value_dim,
                                                   value_offset, space_offset)

        # Update offsets for value- and space dimension
        value_offset += data["value_dim"]
        space_offset += data["space_dim"]

    return code


def _interpolate_vertex_values_element(data, dim, total_value_dim,
                                       value_offset=0, space_offset=0):

    # Extract vertex values for all basis functions
    vertex_values = data["basis_values"]
    value_dim = data["value_dim"]
    space_dim = data["space_dim"]
    mapping = data["mapping"]

    # Create code for each value dimension:
    code = ""
    for k in range(value_dim):

        # Create code for each vertex x_j
        for (j, values_at_vertex) in enumerate(vertex_values):

            if value_dim == 1: values_at_vertex = [values_at_vertex]

            # Map basis values according to element mapping
            if mapping is "affine":
                components = values_at_vertex[k]

            elif mapping == "contravariant piola":
                contraction = lambda i: inner_product([J(k, l) for l in range(dim)],
                                                      [values_at_vertex[l][i] for l in range(dim)])
                components = [multiply([invdetJ, contraction(i)]) for i in range(space_dim)]

            elif mapping == "covariant piola":
                components = [inner_product([Jinv(k, l) for l in range(dim)],
                                            [values_at_vertex[l][i] for l in range(dim)])
                              for i in range(space_dim)]
            else:
                raise Exception, "No such mapping: %s accepted" % mapiping

            # Contract coefficents and basis functions
            dof_values = [component("dof_values", i + space_offset) for i in range(space_dim)]
            value = inner_product(dof_values, components)

            # Assign value to correct vertex
            index = j*total_value_dim + (k + value_offset)
            code += assign(component("vertex_values", index), value)

    return code
