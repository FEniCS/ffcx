"DOLFIN output format."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-14"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

def compile(products, A0s, ranks):
    "Generate code for DOLFIN."
    print "Compiling multi-linear form for C++ (DOLFIN)."
    __write_header()
    __write_geometry_tensor(products, ranks)
    __write_element_tensor(A0s, ranks)
    __write_footer()
    return

def __write_header():
    "Write header for DOLFIN."
    print "class MyPDE : public PDE"
    print "{"
    print "public:"
    print "    void interiorElementMatrix(real** A) const"
    print "    {"
    return

def __write_footer():
    "Write footer for DOLFIN."
    print "    }"
    print "};"
    return

def __write_geometry_tensor(products, ranks):
    "Write expressions for computation of geomety tensor."
    print "        // Compute geometry tensors"
    for i in range(len(ranks)):
        for a in ranks[i].indices1:
            element = "G" + str(i) + "_" + "".join([str(index) for index in a])
            print "        real " + element + " = "
    print
    return

def __write_element_tensor(A0s, ranks):
    """Write expressions for computation of element tensor as
    the product of the geometry tensor and reference tensor."""
    print "        // Compute element tensor"
    for i in ranks[0].indices1: # All primary ranks are equal
        element = "A[" + "][".join([str(index) for index in i]) + "]"
        print "        " + element + " = "
    return
