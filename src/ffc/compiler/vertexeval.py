__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-20 -- 2006-05-30"
__copyright__ = "Copyright (C) 2005-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.dualbasis import *
from FIAT.shapes import *

# FFC modules
from declaration import *

# FIXME: Should not be DOLFIN-specific
format = { ("entity", 2, 0) : lambda i : "cell.connections(0)[%d]" % i,
           ("entity", 2, 1) : lambda i : "cell.connections(1)[%d]" % i,
           ("entity", 2, 2) : lambda i : "cell.index()",
           ("entity", 2, 3) : lambda i : "not defined",
           ("entity", 3, 0) : lambda i : "cell.connections(0)[%d]" % i,
           ("entity", 3, 1) : lambda i : "cell.connections(1)[%d]" % i,
           ("entity", 3, 2) : lambda i : "cell.connections(2)[%d]" % i,
           ("entity", 3, 3) : lambda i : "cell.index()",
           ("num",    2, 0) : "mesh.numVertices()",
           ("num",    2, 1) : "mesh.topology().size(1)",
           ("num",    2, 2) : "mesh.numCells()",
           ("num",    2, 3) : "not defined",
           ("num",    3, 0) : "mesh.numVertices()",
           ("num",    3, 1) : "mesh.topology().size(1)",
           ("num",    3, 2) : "mesh.numFaces()",
           ("num",    3, 3) : "mesh.numCells()",
           ("check",  2, 0) : lambda i : "not defined",
           ("check",  2, 1) : lambda i : "cell.edgeAlignment(%d)" % i,
           ("check",  2, 2) : lambda i : "not defined",
           ("check",  2, 3) : lambda i : "not defined",
           ("check",  3, 0) : lambda i : "not defined",
           ("check",  3, 1) : lambda i : "cell.edgeAlignment(%d)" % i,
           ("check",  3, 2) : lambda i : "cell.faceAlignment(%d)" % i,
           ("check",  3, 3) : lambda i : "not defined" }

def compute_offset(nodes_per_entity, num_dims, dim):
    "Compute increment for global offset."
    increment = None
    if nodes_per_entity[dim] == 1:
        increment = format[("num", num_dims - 1, dim)]
    elif nodes_per_entity[dim] > 1:
        increment = "%d*%s" % (nodes_per_entity[dim], format[("num", num_dims - 1, dim)])
    return increment

def write_offset(increment, component):
    "Write new offset for global nodes."
    if component == 0:
        name = "int offset"
        value = increment
    else:
        name = "offset"
        value = "offset + " + increment
    return Declaration(name, value)

def write_component(component):
    "Write component value."
    if component == 0:
        name = "vertex_nodes[%d]" % component
        value = "vertex"
    else:
        name = "vertex_nodes[%d]" % component
        value = "offset + vertex"
    return Declaration(name, value)

def compute_vertexeval(element, component_offset):
    "Compute node map for given element."

    # Get shape and dual basis from FIAT
    shape = element.fiat_shape
    dual_basis = element.fiat_dual
    
    # Get entity IDs
    entity_ids = dual_basis.entity_ids
        
    # Number of topological dimensions
    num_dims = element.shapedim() + 1

    # Count the nodes associated with each entity
    nodes_per_entity = []
    for dim in range(num_dims):
        if entity_ids[dim]:
            nodes_per_entity += [len(entity_ids[dim][0])]
        else:
            nodes_per_entity += [0]

    # Get the number of vector components
    num_components = dual_basis.num_reps

    # Iterate over vector components
    declarations = []
    for component in range(num_components):

        # Compute global component (handles mixed elements)
        global_component = component_offset + component

        # Write component value
        declarations += [write_component(global_component)]
        
        # Iterate over topological dimensions to compute total offset
        increment = ""
        for dim in range(num_dims):

            # Add to global offset
            if nodes_per_entity[dim] > 0:
                if increment == "":
                    increment = compute_offset(nodes_per_entity, num_dims, dim)
                else:
                    increment += " + " + compute_offset(nodes_per_entity, num_dims, dim)

        # Write offset
        declarations += [write_offset(increment, global_component)]

    # Increase component offset
    component_offset += num_components

    #for declaration in declarations:
    #    print declaration.name + " = " + declaration.value

    return (declarations, component_offset)

class VertexEval:

    """VertexEval evaluates a given component of a function in the
    finite element space at a given vertex"""

    def __init__(self, elements):
        "Create VertexEval."

        # Make sure we have a list of elements
        if not isinstance(elements, list):
            elements = [elements]

        # Iterate over elements (handles mixed elements)
        self.declarations = []
        component_offset = 0
        for element in elements:
            (declarations, component_offset) = compute_vertexeval(element, component_offset)
            self.declarations += declarations

        # Remove last declaration (additional offset)
        self.declarations.pop()
