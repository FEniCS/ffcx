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
from ffc.fem.mixedelement import *

# Python modules
import math
import numpy

def evaluate_basis(element, format):
    """Evaluate element basisfunction at a point, currently only being implemented/tested
       for triangular Lagrange elements. The value of the basisfunction are computed like
       in FIAT as the dot product of the coefficients (computed at compile time) and
       basisvalues which are dependent on the coordinate and thus have to be comuted at
       run time."""

# Supported:
# 2D
# Lagrange                + mixed
# Discontinuous Lagrange  + mixed
# Crouzeix-Raviart        + mixed
# Brezzi-Douglas-Marini   + mixed
# 3D

# Not supported:
# 2D
# Raviart-Thomas ? (not tested since it is broken in FFC, but should work)
# Nedelec (broken?)
#
# 3D
# Lagrange
# Discontinuous Lagrange
# Crouzeix-Raviart
# Raviart-Thomas
# Brezzi-Douglas-Marini (BDM)
# Nedelec (broken?)

# To be fixed:
# check that the code is language-independent
# clean up code

    code = []

    code += [format["comment"]("Not implemented.. yet")]

    # Check if we have just one element
    if (element.num_sub_elements() == 1):
        code += generate_element_code(element, format)

    # If the element is vector valued or mixed
    else:
        code += generate_cases(element, format)

    return code

def generate_element_code(element, format):

    code = []

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

    # Get shape of coefficients
    shape = numpy.shape(coefficients)

    # Scalar valued basis element [Lagrange, Discontinuous Lagrange, Crouzeix-Raviart]
    if (len(shape) == 2):
        num_components = 1
        poly_dim = shape[1]
        coefficients = [coefficients]

    # Vector valued basis element [Raviart-Thomas, Brezzi-Douglas-Marini (BDM)]
    elif (len(shape) == 3):
        num_components = shape[1]
        poly_dim = shape[2]
        coefficients = numpy.transpose(coefficients, [1,0,2])

    # ???
    else:
        raise RuntimeError(), "These coefficients have a strange shape!"

    # Get the number of dofs from element
    num_dofs = element.space_dimension()

    code += [format["comment"]("Table(s) of coefficients")]

    # Generate tables for each component
    for i in range(num_components):

        # Extract coefficients for current component
        coeffs = coefficients[i]

        # Declare varable name for coefficients
        name = format["table declaration"] + "coefficients%d[%d][%d]" %(i, num_dofs, poly_dim,)

        # Generate array of values
        value = "\\\n" + format_block_begin
        rows = []
        for j in range(num_dofs):
            rows += [format_block_begin + ", ".join([format["floating point"](coeffs[j,k])\
                     for k in range(poly_dim)]) + format_block_end]

        value += ",\n".join(rows)
        value += format_block_end

        code += [(name, value)] + [""]

    return code


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

    code += [format["comment"]("Compute values")]

    # Get coefficients from basis functions, computed by FIAT at compile time
    coefficients = element.basis().coeffs

    # Get shape of coefficients
    shape = numpy.shape(coefficients)

    # Scalar valued basis element [Lagrange, Discontinuous Lagrange, Crouzeix-Raviart]
    if (len(shape) == 2):

        # Reset value as it is a pointer
        code += [("*values", "0.0")]

        # Loop dofs to generate dot product, 3D ready
        code += [format["loop"]("j", "j", element.space_dimension(), "j")]
        code += [indent(format["add equal"]("*values","coefficients0[i][j]*basisvalues[j]"),2)]

    # Vector valued basis element [Raviart-Thomas, Brezzi-Douglas-Marini (BDM)]
    elif (len(shape) == 3):
        num_components = shape[1]
        poly_dim = shape[2]

        # Reset value as it is a pointer
        code += [("values[%d]" %(i), "0.0") for i in range(num_components)]

        # Loop dofs to generate dot product, 3D ready
        code += [format["loop"]("j", "j", element.space_dimension(), "j")]
        code += [format["block begin"]]

        code += [indent(format["add equal"]("values[%d]" %(i),\
                 "coefficients%d[i][j]*basisvalues[j]" %(i)),2) for i in range(num_components)]

        code += [format["block end"]]

    # ???
    else:
        raise RuntimeError(), "These coefficients have a strange shape!"

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

def generate_cases(element, format):
    "Generate cases in the event of vector elements or mixed elements"

    # Prefetch formats to speed up code generation
    format_block_begin = format["block begin"]
    format_block_end = format["block end"]

#    print "elements: ", element

    # Extract basis elements
    elements = extract_elements(element)
#    print "extracted elements: ", elements

    code, unique_elements = element_types(elements, format)

#    print "unique_elements: ",unique_elements

    code += dof_map(elements, format)

    num_unique_elements = len(unique_elements)
    if (num_unique_elements > 1):
        code += [format["switch"]("element")]
        code += [format_block_begin]

        for i in range(len(unique_elements)):
            code += [format["case"](i)]
            code += [format_block_begin]
            element = unique_elements[i]
            code += generate_element_code(element, format)
            code += [format["break"]]
            code += [format_block_end]

        code += [format_block_end]

    else:
        element = unique_elements[0]
        code += generate_element_code(element, format)

    return code

def extract_elements(element):
    """This function extracts the individual elements from vector elements and mixed elements.
    Example, the following mixed element:

    element1 = FiniteElement("Lagrange", "triangle", 1)
    element2 = VectorElement("Lagrange", "triangle", 2)

    element  = element2 + element1

    has the structure: mixed-element[mixed-element[Lagrange order 2, Lagrange order 2], Lagrange order 1]

    This function returns the list of basis elements:
    elements = [Lagrange order 2, Lagrange order 2, Lagrange order 1]"""

    elements = [element.sub_element(i) for i in range(element.num_sub_elements())]
    mixed = True
    while (mixed == True):
        mixed = False
        for i in range(len(elements)):
            sub_element = elements[i]
            if isinstance(sub_element, MixedElement):
                mixed = True
                elements.pop(i)
                for j in range(sub_element.num_sub_elements()):
                    elements.insert(i+j, sub_element.sub_element(j))

    return elements

def element_types(elements, format):
    code = []

    # Prefetch formats to speed up code generation
    format_block_begin = format["block begin"]
    format_block_end = format["block end"]

    unique_elements = []
    types = [0]

    for i in range(1, len(elements)):
        unique = True
        element = elements[i]
        elem_type = len(unique_elements)
        for j in range(elem_type):
            if (element.signature() == unique_elements[j].signature()):
                unique = False
                elem_type = j
                break
        if unique:
            unique_elements += [element]
        types += [elem_type]

    # Declare element types and tabulate
    name = format["const uint declaration"] + "element_types[%d]" %(len(elements),)
    value = format_block_begin
    value += ", ".join(["%d" %(element_type) for element_type in types])
    value += format_block_end
    code += [(name, value)]
    return (code, unique_elements)


def dof_map(elements, format):

    code = []

    # Prefetch formats to speed up code generation
    format_block_begin = format["block begin"]
    format_block_end = format["block end"]

    # Declare variable for dof map
    code += [(format["uint declaration"] + "element",0)]
#    code += [(format["uint declaration"] + "dof",0)]
    code += [(format["uint declaration"] + "tmp",0)]

    # Declare dofs_per_element variable and tabulate
    name = format["const uint declaration"] + "dofs_per_element[%d]" %(len(elements),)
    value = format_block_begin
    value += ", ".join(["%d" %(element.space_dimension()) for element in elements])
    value += format_block_end
    code += [(name, value)]

    # Loop elements
    code += [format["loop"]("j", "j", len(elements), "j")]
    code += [format_block_begin]
    # if
    code += ["if (tmp +  dofs_per_element[j] > i)" %()]
    code += [format_block_begin]
    code += [("i", "i - tmp")]
#    code += [("dof", "i - tmp")]
    code += [("element", "element_types[j]")]
    code += [format["break"]]
    code += [format_block_end]
    # else
    code += ["else"]
    code += [format_block_begin]
    code += [format["add equal"]("tmp","dofs_per_element[j]")]
    code += [format_block_end]

    # end loop    
    code += [format_block_end]

    return code


