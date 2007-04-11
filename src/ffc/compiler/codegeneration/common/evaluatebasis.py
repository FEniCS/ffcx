"""Code generation for evaluation of finite element basis values. This module generates
   code which is more or less a C++ representation of FIAT code. More specifically the
   functions from the modules expansion.py and jacobi.py are translated into C++"""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-04-04 -- 2007-04-10"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.constants import *

# FFC fem modules
from ffc.fem.finiteelement import *

# Python modules
import math

def evaluate_basis(element, format):
    """Evaluate element basisfunction at a point, currently only being implemented/tested
       for triangular Lagrange elements. The value of the basisfunction are computed like
       in FIAT as the dot product of the coefficients (computed at compile time) and
       basisvalues which are dependent on the coordinate and thus have to be comuted at
       run time."""

# The current code is working for Lagrange elements of any order
# To be fixed:
# support vector Lagrange 2D
# support Lagrange/vector Lagrange in 3D
# support other element types (also vector valued) in 2D and 3D (hopefully this is working without any additional fixes)
# support mixed elements of all types in 2D and 3D, since vector valued elements are special cases of
#   mixed elements this should be working without to many changes
# check that the code is language-independent

    code = []

    code += [format["comment"]("Not implemented.. yet")]

    # Tabulate coefficients
    code += tabulate_coefficients(element, format)

    # Get coordinates and generate map
    code += generate_map(element, format)

    # Compute scaling 1/2(1-y)^n
    code += compute_scaling(element, format)

    # Compute auxilliary functions
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

    # Prefetch formats to speed up code generation
    format_block_begin  = format["block begin"]
    format_block_end    = format["block end"]

    # Get coefficients from basis functions, computed by FIAT at compile time
    coefficients = element.basis().coeffs

    # Get the number of dofs from element
    num_dofs = element.space_dimension()

    # Declare varable name for coefficients
    code += [format["comment"]("Table of coefficients")]
    name = format["table declaration"] + "coefficients[%d][%d]" %(num_dofs, num_dofs,)

    # Generate array of values
    value = "\\\n" + format_block_begin
    for i in range(num_dofs):
        if i == 0:
            value += format_block_begin
        else:
            value += " " + format_block_begin
        for j in range(num_dofs):
            value += format["floating point"](coefficients[i,j])
            if j == num_dofs - 1:
                value += format_block_end
            else:
                value += ", "
        if i == num_dofs - 1:
            value += format_block_end
        else:
            value += ",\n"

    code += [(name, value)]

    return code + [""]


def generate_map(element, format):
    "Generates map from reference triangle/tetrahedron to reference square/cube"

    code = []

    # Prefetch formats to speed up code generation
    format_comment      = format["comment"]
    format_float        = format["float declaration"]
    format_coordinates  = format["coordinate access"]

    # Code snippets reproduced from FIAT: expansions.py: eta_triangle(xi) & eta_tetrahedron(xi)
    eta_triangle = ["if (%s < %s)" % (format["absolute value"]("y - 1.0"), format["floating point"](FFC_EPSILON),),\
    (indent("x",2), -1.0), "else", (indent("x",2), "2.0 * (1.0 + x)/(1.0 - y) - 1.0")]
    eta_tetrahedron = [format_comment("Mapping from a tetrahedron not implemented yet")]

    # List of coordinate declarations, 3D ready
    coordinates = [(format_float + "x", format_coordinates(0)), \
    (format_float + "y", format_coordinates(1)), (format_float + "z", format_coordinates(2))]

    # Dictionaries, 3D ready
    reference = {2:"square", 3:"cube"}
    mappings = {2:eta_triangle, 3:eta_tetrahedron}

    # Generate code, 3D ready
    code += [format_comment("Get coordinates")]
    code += [coordinates[i] for i in range(element.cell_shape())] + [""]
    code += [format_comment("Map coordinates to the reference %s" % (reference[element.cell_shape()]))]
    code += mappings[element.cell_shape()]

    return code + [""]

def compute_scaling(element, format):
    "Generate the scaling, currently only for 2D"

    # From FIAT: expansions.py: make_scalings(etas)

    code = []

    # Get the element degree
    degree = element.degree()

    # Scale factor, for triangles 1/2*(1-y)^i i being the order of the element
    scale_factor = "(0.5 - 0.5 * y)"

    code += [format["comment"]("Generate scalings")]

    # Declare scaling variable
    name = format["const float declaration"] + "scalings[%d]" %(degree+1,)
    value = format["block begin"] + "1.0"
    if degree > 0:
        value += ", " + ", ".join(["scalings[%d]*%s" %(i-1,scale_factor) for i in range(1, degree+1)])
    value += format["block end"]

    code += [(name, value)]

    return code + [""]

def compute_psitilde_a(element, format):
    "Currently only 2D"

    # From FIAT jacobi.py: jacobi.eval_jacobi(0,0,n,eta1)
    # The format is:
    # psitilde_a[0] = 1.0
    # psitilde_a[1] = a + b * x
    # psitilde_a[n] = a * psitilde_a[n-1] + b * psitilde_a[n-1] * x + c * psitilde_a[n-2]
    # where a, b and c are coefficients computed by eval_jacobi_batch and n is the element degree

    code = []

    # Prefetch formats to speed up code generation
    format_float = format["floating point"]

    # Get the element degree
    degree = element.degree()


    code += [format["comment"]("Compute psitilde_a")]

    # Declare variable
    name = format["table declaration"] + "psitilde_a[%d]" %(degree+1,)

    # Get coefficients
    coeff = eval_jacobi_batch(0,0,degree)
    var = [format_float, lambda i,j: format_float(i) +" * psitilde_a[%d]" %j]
    var1 = [""," * x",""]

    value = format["block begin"] + "1.0"
    if degree > 0:
        if coeff[0][1] < 0:
            value += ", %s" % " ".join([var[0](coeff[0][i])+var1[i] for i in range(2) if coeff[0][i]!=0.0])
        else:
            value += ", %s" % " + ".join([var[0](coeff[0][i])+var1[i] for i in range(2) if coeff[0][i]!=0.0])

        for k in range(2, degree+1):
            signs = [""]
            for j in range(1,3):
                if coeff[k-1][j] < 0:
                    signs += [" - "]
                    coeff[k-1][j] = abs(coeff[k-1][j])
                else:
                    if coeff[k-1][j-1] != 0.0:
                        signs += [" + "]
                        coeff[k-1][j] = abs(coeff[k-1][j])
                    else:
                        signs += [""]
                        coeff[k-1][j] = abs(coeff[k-1][j])

            value += ", %s" % "".join([signs[i] + var[1](coeff[k-1][i], k-1-i*(i-1)*0.5)+var1[i]\
                                       for i in range(3) if coeff[k-1][i] != 0.0])
    value += format["block end"]

    code += [(name, value)]

    return code + [""]

def compute_psitilde_b(element, format):

    # From FIAT jacobi.py: jacobi.eval_jacobi_batch(2*i+1,0,n-i,eta2s) for i in range(0,degree+1)
    # The format is:
    # psitilde_bs_0[0] = 1.0
    # psitilde_bs_0[1] = a + b * y
    # psitilde_bs_0[n] = a * psitilde_bs_0[n-1] + b * psitilde_bs_0[n-1] * x + c * psitilde_bs_0[n-2]
    # psitilde_bs_(n-1)[0] = 1.0
    # psitilde_bs_(n-1)[1] = a + b * y
    # psitilde_bs_n[0] = 1.0
    # where a, b and c are coefficients computed by eval_jacobi_batch
    # n is the element degree + 1

    code = []

    # Prefetch formats to speed up code generation
    format_float = format["floating point"]

    # Get the element degree
    degree = element.degree()


    code += [format["comment"]("Compute psitilde_bs")]
    var = [format_float, lambda i,j,k: format_float(i) +" * psitilde_bs_%d[%d]" %(j,k)]
    var1 = [""," * y",""]

    for i in range(0, degree + 1):
        # Declare variable
        name = format["table declaration"] + "psitilde_bs_%d[%d]" %(i, degree+1-i,)

        # Get coefficients
        coeff = eval_jacobi_batch(2*i+1,0,degree-i)
        value = format["block begin"] + "1.0"
        if degree - i > 0:
            if coeff[0][1] < 0:
                value += ", %s" % " ".join([var[0](coeff[0][j])+var1[j] for j in range(2) if coeff[0][j]!=0.0])
            else:
                value += ", %s" % " + ".join([var[0](coeff[0][j])+var1[j] for j in range(2) if coeff[0][j]!=0.0])

            for k in range(2, degree - i + 1):
                signs = [""]
                for j in range(1,3):
                    if coeff[k-1][j] < 0:
                        signs += [" - "]
                        coeff[k-1][j] = abs(coeff[k-1][j])
                    else:
                        if coeff[k-1][j-1] != 0.0:
                            signs += [" + "]
                            coeff[k-1][j] = abs(coeff[k-1][j])
                        else:
                            signs += [""]
                            coeff[k-1][j] = abs(coeff[k-1][j])

                value += ", %s" % "".join([signs[j] + var[1](coeff[k-1][j], i, k-1-j*(j-1)*0.5)\
                                           + var1[j] for j in range(3) if coeff[k-1][j] != 0.0])
        value += format["block end"]

        code += [(name, value)]

    return code + [""]

def compute_basisvalues(element, format):

    code = []
    code += [format["comment"]("Compute basisvalues")]
    dofs = element.space_dimension()

    # Declare variable
    name = format["table declaration"] + "basisvalues[%d]" %(dofs,)

    value = format["block begin"]
    var = []
    for k in range(0,element.degree()+1):
        for i in range(0,k+1):
            ii = k-i
            jj = i
            factor = math.sqrt( (ii+0.5)*(ii+jj+1.0) )
            var += ["psitilde_a[%d] * scalings[%d] * psitilde_bs_%d[%d] * " %(ii,ii,ii,jj)
                    + format["floating point"](factor)]

    value += ", ".join(var)
    value += format["block end"]

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
    code += [format["loop begin"]("j", "j", element.space_dimension(), "j")]
    code += [indent(format["add equal"]("*values","coefficients[i][j]*basisvalues[j]"),2)]
    code += [format["loop end"]]

    return code

def eval_jacobi_batch(a,b,n):
    """Evaluates all jacobi polynomials with weights a,b
    up to degree n. Returns coefficients in an array wher rows correspond to the Jacobi
    polynomials and the columns correspond to the coefficient."""

    result = []
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
            result += [[a2, a3,-a4]]
    return result










