"DOLFIN output format."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-14"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

EPSILON = 3e-16

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
    for j in range(len(ranks)):
        for a in ranks[j].indices1:
            name_g = __name_g(j, a)
            value_g = __value_g(products[j], a, ranks[j].indices2)
            print "        real " + name_g + " = " + value_g + ";"
    print
    return

def __write_element_tensor(A0s, ranks):
    """Write expressions for computation of element tensor as
    the product of the geometry tensor and reference tensor."""
    print "        // Compute element tensor"
    for i in ranks[0].indices0: # All primary ranks are equal
        name_a = __name_a(i)
        value_a = __value_a(A0s, ranks, i)
        print "        " + name_a + " = " + value_a + ";"
    return

def __name_g(j, a):
    "Return name of current element of current geometry tensor."
    return "G" + str(j) + "_" + "".join([str(index) for index in a])

def __value_g(product, a, indices2):
    "Return value of current element of current geometry tensor."
    value = ""
    # Include constant
    if product.constant == -1.0:
        value = "-"
    elif not product.constant == 1.0:
        value = str(product.constant) + "*"
    # Expression handled differently if we have a sum
    if indices2:
        if not product.constant == 1.0:
            value += "("
        terms = ["*".join([__value_transform(t, a, b) \
                          for t in product.transforms]) for b in indices2]
        value += " + ".join(terms)
        if not product.constant == 1.0:
            value += ")"
    else:
        value += "*".join([__value_transform(t, a, []) for t in product.transforms])
    # FIXME: Include coefficients
    return value

def __name_a(i):
    "Return name of current element of element tensor."
    return "A[" + "][".join([str(index) for index in i]) + "]"

def __value_a(A0s, ranks, i):
    """Return value of current element of element tensor. Sum over
    all products and for each product compute the tensor product."""
    value = ""
    for j in range(len(A0s)):
        if ranks[j].indices1:
            for a in ranks[j].indices1:
                element = A0s[j][i + a]
                if abs(element) > EPSILON:
                    if value and element < 0.0:
                        value += " - " + str(-element) + "*" + __name_g(j, a)
                    elif value:
                        value += " + " + str(element) + "*" + __name_g(j, a)
                    else:
                        value += str(element) + "*" + __name_g(j, a)
        else:
            element = A0s[j][i]
            if value and element < 0.0:
                value += " - " + str(-element)
            elif value:
                value += " + " + str(element)
            else:
                value += str(element)
    if '+' in value or '-' in value:
        return "det*(" + value + ")"
    else:
        return "det*" + value

def __value_transform(transform, a, b):
    "Return value of current transform."
    value = "g"
    value += str(transform.index0([], a, b))
    value += str(transform.index1([], a, b))
    return value
