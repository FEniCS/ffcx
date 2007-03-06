"Code generation for dof map"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-24 -- 2007-02-06"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC fem modules
from ffc.fem.finiteelement import *

def generate_dof_map(dof_map, format):
    """Generate dictionary of code for the given dof map according to
    the given format"""

    code = {}

    # Generate code for signature
    code["signature"] = dof_map.signature()

    # Generate code for needs_mesh_entities
    code["needs_mesh_entities"] = __generate_needs_mesh_entities(dof_map, format)

    # Generate code for global_dimension
    code["global_dimension"] = __generate_global_dimension(dof_map, format)

    # Generate code for local dimension
    code["local_dimension"] = "%d" % dof_map.local_dimension() 

    # Generate code for tabulate_dofs
    code["tabulate_dofs"] = __generate_tabulate_dofs(dof_map, format)

    return code

def __generate_needs_mesh_entities(dof_map, format):
    "Generate code for needs_mesh_entities"

    # Get the number of dofs per dimension
    dofs_per_dimension = dof_map.dofs_per_dimension()

    # Entities needed if at least one dof is associated
    code = [format["bool"](num_dofs > 0) for num_dofs in dofs_per_dimension]

    return code

def __generate_global_dimension(dof_map, format):
    "Generate code for global dimension"

    # Get the number of dofs per dimension
    dofs_per_dimension = dof_map.dofs_per_dimension()

    # Sum the number of dofs for each dimension
    terms = []
    for dim in range(len(dofs_per_dimension)):
        n = dofs_per_dimension[dim]
        if n == 1:
            terms += [format["num entities"](dim)]
        elif n > 1:
            terms += [format["multiply"]([str(n), format["num entities"](dim)])]

    # Special case, no terms
    if len(terms) == 0:
        code = "0"
    else:
        code = format["add"](terms)

    return code

def __generate_tabulate_dofs(dof_map, format):
    "Generate code for tabulate_dofs"

    # Generate code as a list of declarations
    code = []

    # Get entity dofs and dofs per dimension
    entity_dofs = dof_map.entity_dofs()
    dofs_per_dimension = dof_map.dofs_per_dimension()

    #for dim in entity_dofs:
    #    print "dim = " + str(dim) + ": " + str(entity_dofs[dim])

    # Iterate over dimensions
    offset_declared = False
    offset_code = []
    for dim in entity_dofs:

        # Skip dimension if there are no dofs
        if dofs_per_dimension[dim] == 0:
            continue

        # Write offset code
        code += offset_code

        # Iterate over entities in dimension
        for entity in entity_dofs[dim]:

            # Iterate over dofs on entity
            for pos in range(len(entity_dofs[dim][entity])):

                # Get number of current dof
                dof = entity_dofs[dim][entity][pos]

                # Assign dof
                name = format["dofs"](dof)
                value = format["entity index"](dim, entity)

                # Add position on entity if any
                if pos > 0:
                    value = format["add"]([value, "%d" % pos])

                # Add offset if any
                if offset_declared:
                    value = format["add"]([format["offset access"], value])

                # Add declaration
                code += [(name, value)]

        # Add to offset if needed
        if dofs_per_dimension[dim] > 0:
            if dofs_per_dimension[dim] > 1:
                value = format["multiply"](["%d" % dofs_per_dimension[dim], format["num entities"](dim)])
            else:
                value = format["num entities"](dim)
            if not offset_declared:
                name = format["offset declaration"]
                offset_declared = True
            else:
                name = format["offset access"]
                value = format["add"]([name, value])
                
            offset_code += [(name, value)]

    return code
