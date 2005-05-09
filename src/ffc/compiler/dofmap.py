__author__ = "Robert C. Kirby (kirby@cs.uchicago.edu) and Anders Logg (logg@tti-c.org)"
__date__ = "2005-05-03 -- 2005-05-09"
__copyright__ = "Copyright (c) 2005 Kirby/Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.dualbasis import *
from FIAT.shapes import *

# FFC modules
from declaration import *

format = { ("entity", 0) : lambda i : "cell.nodeID(%d)" % i,
           ("entity", 1) : lambda i : "cell.edgeID(%d)" % i,
           ("entity", 2) : lambda i : "cell.id()",
           ("entity", 3) : lambda i : "cell.id()",
           ("num", 0)    : "mesh.noNodes()",
           ("num", 1)    : "mesh.noEdges()",
           ("num", 2)    : "mesh.noFaces()",
           ("num", 3)    : "mesh.noCells()",
           ("check", 0)  : lambda i : "",
           ("check", 1)  : lambda i : "cell.edgeAligned(%d)" % i,
           ("check", 2)  : lambda i : "",
           ("check", 3)  : lambda i : ""  }

class DofMap:

    """A DofMap maps the dofs (degrees of freedom) on a local element
    to global degrees of freedom."""

    def __init__(self, shape, dualbasis):
        "Create DofMap."

        print "---------------------------------------------------"
        print "Creating dof map (experimental)"
        print ""

        print "entity ids:   " + str(dualbasis.entity_ids)

        # Number of topological dimensions
        num_dims = dimension(shape) + 1

        # Count the entities associated with each topological dimension
        num_entities = [len(entity_range(shape, dim)) for dim in range(num_dims)]
        print "num entities: " + str(num_entities)

        # Count the nodes associated with each entity
        num_nodes = [len(dualbasis.getNodeIDs(dim)[0]) for dim in range(num_dims)]
        print "num nodes:    " + str(num_nodes)

        print ""
        
        self.declarations = []
        current = 0
        offset = []
        
        # Iterate over topological dimensions
        for dim in range(num_dims):
            # Iterate over entities for current dimension
            for entity in range(num_entities[dim]):
                if num_nodes[dim] == 1:
                    # Map the single node associated with the current entity
                    name = "dofs[%d]" % current
                    value = " + ".join(offset + [format[("entity", dim)](entity)])
                    self.declarations += [Declaration(name, value)]
                    current += 1
                elif num_nodes[dim] > 1:
                    # Iterate over the nodes associated with the current entity
                    start = "%d*%s" % (num_nodes[dim], format[("entity", dim)](entity))
                    for node in range(num_nodes[dim]):
                        if dim < (num_dims - 1):
                            local = "( %s ? %d : %d )" % (format[("check", dim)](entity), node, num_nodes[dim] - 1 - node)
                        else:
                            local = "%d" % node
                        name = "dofs[%d]" % current
                        value = " + ".join(offset + [start] + [local])
                        self.declarations += [Declaration(name, value)]
                        current += 1
            # Add to offset
            if num_nodes[dim] == 1:
                offset += [format[("num", dim)]]
            elif num_nodes[dim] > 1:
                offset += ["%d*%s" % (num_nodes[dim], format[("num", dim)])]

        for declaration in self.declarations:
            print declaration.name + " = " + declaration.value

        print "---------------------------------------------------"
