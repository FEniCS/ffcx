"DOLFIN output format."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-14"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

def compile(A0):
    "Generate code for DOLFIN."
    print "Compiling multi-linear form for C++ (DOLFIN)."
    __write_header()
    __write_geometry_tensor()
    __write_element_tensor(A0)
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

def __write_geometry_tensor():
    "Write expressions for computation of geomety tensor."
    print "        G00 = "
    return

def __write_element_tensor(A0):
    """Write expressions for computation of element tensor as
    the product of the geometry tensor and reference tensor."""
    print "        A[0][0] = "
    return
