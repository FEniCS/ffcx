"""Code generation for finite element basis values. This module generates code which is more
   or less a C++ representation of FIAT code. More specifically the functions from
   the modules expansion.py and jacobi.py are translated into C++"""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-04-04 -- 2007-04-04"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# FFC fem modules
from ffc.fem.finiteelement import *

def evaluate_basis(element, format):
    """Evaluate element basisfunction at a point, currently only being implemented/tested
       for triangular Lagrange elements. The value of the basisfunction are computed like
       in FIAT as the dot product of the coefficients (computed at compile time) and
       basisvalues which are dependent on the coordinate and thus have to be comuted at
       run time."""

    code = []

    code += [format["comment"]("Not implemented.. yet")]

    # Tabulate coefficients
    code += tabulate_coefficients(element, format)

    # Get coordinates and generate map
    code += generate_map(element, format)

    # Compute scaling 1/2(1-y)^n
    code += compute_scaling(element, format)

    # Don't know exactly what these functions do/are
    code += compute_psitilde_a(element, format)
    code += compute_psitilde_b(element, format)

    # Compute the basisvalues
    code += compute_basisvalues(element, format)

    # Compute the value of the basisfunction as the dot product of the coefficients
    # and basisvalues
    code += dot_product(element, format)

    return code

def tabulate_coefficients(element, format):

    code = []

    # Get coefficients from basis functions, computed by FIAT at compile time
    coefficients = element.basis().coeffs

    # Get the number of dofs from element
    num_dofs = element.space_dimension()

    # Declare varable name for coefficients
    code += [format["comment"]("Table of coefficients")]
    name = "const static double coefficients[%d][%d]" %(num_dofs, num_dofs,)

    # Generate array of values
    value = "\\\n{"
    for i in range(num_dofs):
        if i == 0:
            value += "{"
        else:
            value += " {"
        for j in range(num_dofs):
            value += "%f" %(coefficients[i,j],)
            if j == num_dofs - 1:
                value += "}"
            else:
                value += ", "
        if i == num_dofs - 1:
            value += "}"
        else:
            value += ",\n"

    code += [(name, value)]

    return code + [""]


def generate_map(element, format):
    "Generates map from reference triangle/tetrahedron to reference square/cube"

    code = []

    # Code snippets reproduced from FIAT: expansions.py: eta_triangle(xi) & eta_tetrahedron(xi)
    eta_triangle = ["if ((y - 1.0) < DOLFIN_EPS)", indent("x = -1.0;",2), "else",\
                    indent("x = 2.0 * (1.0 + x)/(1.0 - y);",2)]
    eta_tetrahedron = [format["comment"]("Mapping from a tetrahedron not implemented yet")]

    # List of coordinate declarations, 3D ready
    coordinates = ["double x = coordinates[0];", "double y = coordinates[1];", "double z = coordinates[2];"]

    # Dictionaries, 3D ready
    reference = {2:"square", 3:"cube"}
    mappings = {2:eta_triangle, 3:eta_tetrahedron}

    # Generate code, 3D ready
    code += [format["comment"]("Get coordinates")]
    code += [coordinates[i] for i in range(element.cell_shape())] + [""]
    code += [format["comment"]("Map coordinates to the reference %s" % (reference[element.cell_shape()]))]
    code += mappings[element.cell_shape()]

    return code + [""]

def compute_scaling(element, format):
    "Generate the scaling, currently only for 2D"

    code = []

    # Get the element degree
    degree = element.degree()

    # From FIAT: expansions.py: make_scalings(etas)

    # Scale factor, for triangles 1/2*(1-y)^i i being the order of the element
    scale_factor = "(0.5 - 0.5 * y)"

    code += [format["comment"]("Generate scaling")]

    # Declare scaling variable
    name = "const double scaling[%d]" %(degree+1,)
    value = "{1.0"
    if degree > 0:
        value += ", " + ", ".join(["scaling[%d]*%s" %(i-1,scale_factor) for i in range(1, degree+1)])
    value += "}"

    code += [(name, value)]

    return code + [""]

def compute_psitilde_a(element, format):
    "Currently only 2D"

    code = []

    # Get the element degree
    degree = element.degree()

    # From FIAT jacobi.py: jacobi.eval_jacobi(0,0,n,eta1)

    code += [format["comment"]("Compute psitilde_a")]

    # Declare variable
    name = "const double psitilde_a[%d]" %(degree+1,)

    # Get coefficients
    coeff = eval_jacobi_batch(0,0,degree)
    var = [["%f", "%f*x"], ["(%f * x) * psitilde_a[%d]","%f * psitilde_a[%d]"]]

    value = "{1.0"
    if degree > 0:
        value += ", %s" % " + ".join([var[0][i] % (coeff[0][i],) for i in range(2) if coeff[0][i]!=0.0])
        for k in range(2, degree+1):
            if coeff[k-1][1] < 0:
                value += ", %s" %" ".join([var[1][i] % (coeff[k-1][i],k-1-i,) \
                         for i in range(2) if coeff[k-1][i]!=0.0])
            else:
                value += ", %s" %" + ".join([var[1][i] % (coeff[k-1][i],k-1-i) \
                         for i in range(2) if coeff[k-1][i]!=0.0])
    value += "}"

    code += [(name, value)]

    return code + [""]

def compute_psitilde_b(element, format):

    code = []

    # Get the element degree
    degree = element.degree()

    # From FIAT jacobi.py: jacobi.eval_jacobi_batch(2*i+1,0,n-i,eta2s) for i in range(0,degree+1)

#    psitilde_bs = [ jacobi.eval_jacobi_batch(2*i+1,0,n-i,eta2s) \
#		    for i in range(0,n+1) ]

    code += [format["comment"]("Compute psitilde_bs")]

    for i in range(0, degree + 1):
        # Declare variable
        name = "const double psitilde_b_%d[%d]" %(i, degree+1-i,)
        value = ""
        code += [(name, value)]

    return code + [""]

def compute_basisvalues(element, format):

    code = []
    code += [format["comment"]("Compute basisvalues")]
    # Declare variable
    name = "const double basisvalues[%d]" %(element.space_dimension(),)
    value = ""

    code += [(name, value)]

    return code + [""]


def dot_product(element, format):
    """This function computes the value of the basisfunction as the dot product of the
       coefficients and basisvalues """

    code = []

    code += [format["comment"]("Compute value")]
    # Reset value as it is a pointer
    code += [("*values", "0.0")]
    # Loop dofs to generate dot product, 3D ready
    code += ["for (unsigned int j = 0; j < %d; j++)" % (element.space_dimension(),)]
    code += [indent(format["add equal"]("*values","coefficients[i][j]*coefficients[i][j];"),2)]
#    code += [indent(format["add equal"]("values","coefficients[i][j]*bvals[j];"),2)]

    return code

def eval_jacobi_batch(a,b,n):
    """Evaluates all jacobi polynomials with weights a,b
    up to degree n. Returns coefficients in an array wher rows correspond to the Jacobi
    polynomials and the columns correspond to the coefficient."""

#    result = numpy.zeros( (n+1, len(xs)),"d" )
#    result = []#[[1.0]]
    if n > 0:
        result = [[0.5 * (a - b), 0.5 * ( a + b + 2.0 )]]
        apb = a + b
        for k in range(2,n+1):
            a1 = 2.0 * k * ( k + apb ) * ( 2.0 * k + apb - 2.0 )
            a2 = ( 2.0 * k + apb - 1.0 ) * ( a * a - b * b )
            a3 = ( 2.0 * k + apb - 2.0 ) * ( 2.0 * k + apb - 1.0 ) * ( 2.0 * k + apb )
            a4 = 2.0 * ( k + a - 1.0 ) * ( k + b - 1.0 ) * ( 2.0 * k + apb )
            a2 = a2 / a1
            a3 = a3 / a1
            a4 = a4 / a1
            result += [[( a2 + a3 ),-a4]]
    return result










