"Code generation for dof map"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-24 -- 2007-04-10"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.utils import *

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

    # Generate code for local_dimension
    code["local_dimension"] = "%d" % dof_map.local_dimension() 

    # Generate code for num_facet_dofs
    code["num_facet_dofs"] = "%d" % dof_map.num_facet_dofs()

    # Generate code for tabulate_dofs
    code["tabulate_dofs"] = __generate_tabulate_dofs(dof_map, format)

    # Generate code for tabulate_facet_dofs
    code["tabulate_facet_dofs"] = __generate_tabulate_facet_dofs(dof_map, format)

    return code

def __generate_needs_mesh_entities(dof_map, format):
    "Generate code for needs_mesh_entities"

    # Get total number of dofs per dimension
    num_dofs_per_dim = dof_map.num_dofs_per_dim()

    # Entities needed if at least one dof is associated
    code = [format["bool"](num_dofs_per_dim[dim] > 0) for dim in range(len(num_dofs_per_dim))]

    return code

def __generate_global_dimension(dof_map, format):
    "Generate code for global dimension"

    # Get total number of dofs per dimension
    num_dofs_per_dim = dof_map.num_dofs_per_dim()
    
    # Sum the number of dofs for each dimension
    terms = []
    for dim in range(len(num_dofs_per_dim)):
        n = num_dofs_per_dim[dim]
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

    # Iterate over sub dofs
    offset_declared = False
    offset_code = []
    local_offset = 0
    for sub_dof_map in range(len(dof_map.entity_dofs())):

        # Get entity dofs for sub dof map
        sub_entity_dofs = dof_map.entity_dofs()[sub_dof_map]

        # Get the number of dofs per dimension for sub dof map
        num_dofs_per_dim = dof_map.num_dofs_per_dim(sub_dof_map)

        # Iterate over dimensions
        num_dofs = 0
        for dim in sub_entity_dofs:

            # Skip dimension if there are no dofs
            if num_dofs_per_dim[dim] == 0:
                continue

            # Write offset code
            code += offset_code

            # Iterate over entities in dimension
            for entity in sub_entity_dofs[dim]:

                # Iterate over dofs on entity
                for pos in range(len(sub_entity_dofs[dim][entity])):

                    # Get number of current dof
                    dof = sub_entity_dofs[dim][entity][pos]

                    # Assign dof
                    name = format["dofs"](local_offset + dof)
                    if num_dofs_per_dim[dim] > 1:
                        value = format["multiply"](["%d" % num_dofs_per_dim[dim], format["entity index"](dim, entity)])
                    else:
                        value = format["entity index"](dim, entity)
                    
                    # Add position on entity if any
                    if pos > 0:
                        value = format["add"]([value, "%d" % pos])

                    # Add offset if any
                    if offset_declared:
                        value = format["add"]([format["offset access"], value])

                    # Add declaration
                    code += [(name, value)]

                    # Count the number of dofs for sub dof map
                    num_dofs += 1

            # Update offset
            if num_dofs_per_dim[dim] > 0:

                # Compute additional offset
                if num_dofs_per_dim[dim] > 1:
                    value = format["multiply"](["%d" % num_dofs_per_dim[dim], format["num entities"](dim)])
                else:
                    value = format["num entities"](dim)

                # Add to previous offset
                if not offset_declared:
                    name = format["offset declaration"]
                    offset_declared = True
                else:
                    name = format["offset access"]
                    value = format["add"]([name, value])
                
                offset_code = [(name, value)]

        # Add to local offset
        local_offset += num_dofs 

    return code

def __generate_tabulate_facet_dofs(dof_map, format):
    "Generate code for tabulate_dofs"

    return ["// Not implemented"]

    # Get the number of facets
    num_facets = dof_map.element().num_facets()
    print num_facets

    # Get incidence
    incidence = dof_map.incidence()

    # Count the number of facets
    incident_entities = num_facets*[{}]

    print ""

    D = 2
    
    # Find out which entities are incident with each facet
    for facet in range(num_facets):
        incident_entities[facet] = [pair[1] for pair in incidence if incidence[pair] and pair[0] == (D - 1, facet)]
        print incident_entities[facet]

    print incident_entities
    

