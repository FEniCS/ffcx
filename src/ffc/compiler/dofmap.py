__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-05-03"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.dualbasis import *
from FIAT.shapes import *

format = { ("entity", 0) : lambda i : "cell.nodeID(%s)" % i,
           ("entity", 1) : lambda i : "cell.edgeID(%s)" % i,
           ("entity", 2) : lambda i : "cell.faceID(%s)" % i,
           ("entity", 3) : lambda i : "cell.id()",
           ("num", 0)    : "mesh.noNodes()",
           ("num", 1)    : "mesh.noEdges()",
           ("num", 2)    : "mesh.noFaces()",
           ("num", 3)    : "mesh.noCells()" }

class DofMap:

    """A DofMap maps the dofs (degrees of freedom) on a local element
    to global degrees of freedom."""

    def __init__(self, shape, dualbasis):
        "Create DofMap."

        print "Creating dof map for " + str(dualbasis)
        print ""

        # Number of topological dimensions
        num_dims = dimension(shape) + 1

        # Count the entities associated with each topological dimension
        num_entities = [len(entity_range(shape, dim)) for dim in range(num_dims)]
        print "num entities:   " + str(num_entities)

        # Count the dofs associated with each topological dimension
        num_dofs = [len(dualbasis.getNodeIDs(dim)[0]) for dim in range(num_dims)]
        print "num dofs    :   " + str(num_dofs)

        # Compute global offsets for each topological dimension
        local_offsets = [0]
        global_offsets = [None]
        for dim in range(1, num_dims):
            # Compute increment
            local_offset = num_dofs[dim - 1] * num_entities[dim - 1]
            global_offset = str(num_dofs[dim - 1]) + "*" + format[("num", dim - 1)]
            # Add previous offset
            local_offset += local_offsets[dim - 1]
            if global_offsets[dim - 1]:
                global_offset = global_offsets[dim - 1] + " + " + global_offset
            # Add to list
            local_offsets += [local_offset]
            global_offsets += [global_offset]
        print "local offsets:  " + str(local_offsets)
        print "global offsets: " + str(global_offsets)
        print ""

        # Generate code for mapping
        self.output = ""
        sum_dofs = 0
        for dim in range(num_dims):
            sum_dofs += num_dofs[dim] * num_entities[dim]
            if dim == 0:
                self.output += "if ( i < %d ) return %s;\n" % \
                               (sum_dofs,
                                format[("entity", dim)]("i - " + str(local_offsets[dim])))
            else:
                self.output += "else if ( i < %d ) return %s + %s;\n" % \
                               (sum_dofs, global_offsets[dim],
                                format[("entity", dim)]("i - " + str(local_offsets[dim])))
        print self.output
                   
        dofs_0 = dualbasis.getNodeIDs(0)
        dofs_1 = dualbasis.getNodeIDs(1)
        dofs_2 = dualbasis.getNodeIDs(2)
        #dofs_3 = dualbasis.getNodeIDs(3)
 
        print "Vertices: " + str(dofs_0)
        print "Edges:    " + str(dofs_1)
        print "Faces:    " + str(dofs_2)
        #print "Cells:    " + str(dofs_3)
