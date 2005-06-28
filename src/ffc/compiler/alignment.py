__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-06-27"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from math import sqrt

def edge_reordering(num_nodes):
    """Compute and return map for reordering of degrees of freedom on
    edges. The return value is a two-dimensional list, indexed by
    (alignment, node number), which maps degrees of freedom from the
    numbering local to the edge of a cell (triangle or tetrahedron) to
    a common local numbering for two cells sharing a common edge. The
    reordering is determined by the alignment of the common edge with
    the cell."""

    return [range(num_nodes), range(num_nodes-1,-1,-1)]

def face_reordering(num_nodes):
    """Compute and return map for reordering of degrees of freedom on
    faces of tetrahedra. The return value is a two-dimensional list,
    indexed by (alignment, node number), which maps degrees of freedom
    from the numbering local to the face of a tetrahedron to a common
    local numbering for two tetrahedra sharing a common face. The
    reordering is determined by the alignment of the common face with
    the tetrahedron."""

    # Compute number of nodes in each direction
    n = int((sqrt(8*num_nodes+1)-3)/2 + 1)
    if not n*(n+1)/2 == num_nodes:
        raise RuntimeError, "Inconsistent number of nodes on face."
    
    # Offsets for rows (n nodes in first row, n-1 in second and so on)
    offsets = [0]
    for i in range(1, n):
        offsets += [offsets[-1] + n + 1 - i]
    print "offsets = " + str(offsets)

    # Now do a nested loop over rows (i) and columns (j) in different
    # directions depending on the alignment
    nodes0 = [0 for i in range(num_nodes)]
    nodes1 = [0 for i in range(num_nodes)]
    nodes2 = [0 for i in range(num_nodes)]
    nodes3 = [0 for i in range(num_nodes)]
    nodes4 = [0 for i in range(num_nodes)]
    nodes5 = [0 for i in range(num_nodes)]

    # alignment == 0
    node = 0
    for i in range(n):
        for j in range(n-i):
            nodes0[offsets[i] + j] = node
            node += 1
    # alignment == 1
    node = 0
    for j in range(n):
        for i in range(n-j):
            nodes1[offsets[i] + j] = node
            node += 1
    # alignment == 2
    node = 0
    for j in range(n):
        for i in range(n-1-j,-1,-1):
            nodes2[offsets[i] + j] = node
            node += 1
    # alignment == 3
    node = 0
    for i in range(n):
        for j in range(n-1-i,-1,-1):
            nodes3[offsets[i] + j] = node
            node += 1
    # alignment == 4
    node = 0
    for j in range(n-1,-1,-1):
        for i in range(j + 1):
            nodes4[offsets[i] + j-i] = node
            node += 1
    # alignment == 5
    node = 0
    for j in range(n-1,-1,-1):
        for i in range(j,-1,-1):
            nodes5[offsets[i] + j-i] = node
            node += 1
            
    print nodes0
    print nodes1
    print nodes2
    print nodes3
    print nodes4
    print nodes5
        
    #for alignment in range(6):
    #    print alignment

    return [nodes0, nodes1, nodes2, nodes3, nodes4, nodes5]
