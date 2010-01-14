"Code generation for interpolate_vertex_values."

__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2009"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-14

from ffc.codesnippets import jacobian
from ffc.cpp import format, remove_unused

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

def interpolate_vertex_values(ir):
    "Generate code for interpolate_vertex_values."

    # Add code for Jacobian if necessary
    code = []
    dim = ir["cell_dim"]
    if ir["needs_jacobian"]:
        code.append(jacobian[dim] % {"restriction": ""})

    # Compute total value dimension for (mixed) element
    total_dim = sum(data["value_dim"] for data in ir["element_data"])

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
        value_offset += data["value_dim"]
        space_offset += data["space_dim"]

    # Remove unused variables. (Not tracking set of used variables in
    # order to keep this code clean. Since generated code is of
    # limited size, this should be ok.)
    clean_code = remove_unused("\n".join(code))
    return clean_code


def _interpolate_vertex_values_element(data, dim, total_value_dim,
                                       value_offset=0, space_offset=0):

    # Extract vertex values for all basis functions
    vertex_values = data["basis_values"]
    value_dim = data["value_dim"]
    space_dim = data["space_dim"]
    mapping = data["mapping"]

    # Map basis values according to element mapping. Assumes single
    # mapping for each (non-mixed) element
    change_of_variables = _change_variables(data["mapping"], dim, space_dim)

    # Create code for each value dimension:
    code = []
    for k in range(value_dim):

        # Create code for each vertex x_j
        for (j, values_at_vertex) in enumerate(vertex_values):

            if value_dim == 1: values_at_vertex = [values_at_vertex]

            # Map basis functions using appropriate mapping
            components = change_of_variables(values_at_vertex, k)

            # Contract coefficients and basis functions
            dof_values = [component("dof_values", i + space_offset)
                          for i in range(space_dim)]
            value = inner_product(dof_values, components)

            # Assign value to correct vertex
            index = j*total_value_dim + (k + value_offset)
            code.append(assign(component("vertex_values", index), value))

    return "\n".join(code)


def _change_variables(mapping, dim, space_dim):

    if mapping is "affine":
        change_of_variables = lambda v, k: v[k]
    elif mapping == "contravariant piola":
        change_of_variables = lambda v, k: [multiply([invdetJ, inner_product([J(k, l) for l in range(dim)],
                                                                             [v[l][i] for l in range(dim)])])
                                            for i in range(space_dim)]
    elif mapping == "covariant piola":
        change_of_variables = lambda v, k: [inner_product([Jinv(k, l) for l in range(dim)],
                                                          [v[l][i] for l in range(dim)])
                                            for i in range(space_dim)]
    else:
        raise Exception, "No such mapping: %s accepted" % mapiping
    return change_of_variables
