"""Code generation for evaluation of derivatives of finite element basis functions. This
module generates code which is more or less a C++ representation of FIAT code. More
specifically the functions from the modules expansion.py and jacobi.py are translated into C++"""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-04-16 -- 2007-04-16"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.constants import *

# FFC fem modules
from ffc.fem.finiteelement import *
from ffc.fem.mixedelement import *

# FFC code generation common modules
import evaluatebasis

# Python modules
import math
import numpy

def evaluate_basis_derivatives(element, format):
    """Evaluate the derivatives of an element basisfunction at a point. The values are
    computed as in FIAT as the dot product of the coefficients (computed at compile time)
    and basisvalues which are dependent on the coordinate and thus have to be computed at
    run time.

    Currently the following elements are supported in 2D and 3D:

    Not supported in 2D or 3D:

    Lagrange                + mixed/vector valued
    Discontinuous Lagrange  + mixed/vector valued
    Crouzeix-Raviart        + mixed/vector valued
    Brezzi-Douglas-Marini   + mixed/vector valued
    Raviart-Thomas ? (not tested since it is broken in FFC, but should work)
    Nedelec (broken?)"""

# To be fixed:
# doc strings
# clean code
# check that the code is language-independent

    code = []

    Indent = evaluatebasis.IndentControl()

    # Get coordinates and generate map
    code += evaluatebasis.generate_map(element, Indent, format)

    # Check if we have just one element
    if (element.num_sub_elements() == 1):
        code += evaluatebasis.dof_map(0, Indent, format)
        code += generate_element_code(element, 0, False, Indent, format)

    # If the element is vector valued or mixed
    else:
        code += generate_cases(element, Indent, format)

    return code

def generate_element_code(element, value_num, vector, Indent, format):
    "Generate code for each basis element"

    code = []

    # Compute basisvalues, from evaluatebasis.py
    code += compute_basisvalues(element, Indent, format)

    # Tabulate coefficients
    code += evaluatebasis.tabulate_coefficients(element, Indent, format)

    code += [Indent.indent(format["comment"]("Interesting part"))]

    # Tabulate coefficients for derivatives
    code += tabulate_dmats(element, Indent, format)

    # Compute the value of the derivatives of the basisfunctions as the dot product of the new coefficients
    # and basisvalues
    code += compute_values(element, Indent, format)

    return code

def generate_cases(element, Indent, format):
    "Generate cases in the event of vector valued elements or mixed elements"

    code = []

    # Prefetch formats to speed up code generation
    format_block_begin = format["block begin"]
    format_block_end = format["block end"]

    # Extract basis elements, and determine number of elements
    elements = evaluatebasis.extract_elements(element)
    num_elements = len(elements)

    # Loop all elements
    code += [Indent.indent(format["loop"]("element", num_elements))]

    code += [Indent.indent(format["comment"]("Switch for each of the basis elements"))]
    code += [Indent.indent(format_block_begin)]
    # Increase indentation
    Indent.increase()

    # Generate switch
    code += [Indent.indent(format["switch"]("element"))]
    code += [Indent.indent(format_block_begin)]
    # Increase indentation
    Indent.increase()

    sum_value_num = 0
    sum_space_dim = 0

    # Generate cases
    for i in range(num_elements):
        code += [Indent.indent(format["case"](i))]
        code += [Indent.indent(format_block_begin)]
        # Increase indentation
        Indent.increase()

        # Get sub element
        basis_element = elements[i]

        # FIXME: This must most likely change for tensor valued elements
        value_dimension = basis_element.value_dimension(0)
        value_num = basis_element.value_dimension(0)
        space_dim = basis_element.space_dimension()

        # Determine if the element has a value, for the given dof
# Not languge-independent
        code += [Indent.indent("if (%d <= i and i <= %d)\n{" % (sum_space_dim, sum_space_dim + space_dim -1))]
        # Increase indentation
        Indent.increase()

        # Generate map from global to local dof
        code += [Indent.indent(format["comment"]("Compute local degree of freedom"))]
        code += evaluatebasis.dof_map(sum_space_dim, Indent, format)

        # Generate code for basis element
        code += generate_element_code(basis_element, sum_value_num, True, Indent, format)

        # Decrease indentation, finish block - end element code
        Indent.decrease()
        code += [Indent.indent(format_block_end)]

        # If the element does not have a value for the given dof, return 0.0
# Not languge-independent
        code += [Indent.indent("else\n{")]
        # Increase indentation
        Indent.increase()
        # Reset values
        code += evaluatebasis.reset_values(value_dimension, sum_value_num, True, Indent, format)
        # Decrease indentation
        Indent.decrease()
        code += [Indent.indent(format_block_end)]

        # End case
        code += [Indent.indent(format["break"])]
        # Decrease indentation
        Indent.decrease()
        code += [Indent.indent(format_block_end)]

        # Increase sum of value dimension, and space dimension
        sum_value_num += value_num
        sum_space_dim += space_dim

    # Decrease indentation, end switch
    Indent.decrease()
    code += [Indent.indent(format_block_end)]

    # Decrease indentation, end loop elements
    Indent.decrease()
    code += [Indent.indent(format_block_end)]

    return code

def compute_basisvalues(element, Indent, format):
    "Code generation from evaluatebasis.py"

    code = []

    # Compute scaling of y and z 1/2(1-y)^n and 1/2(1-z)^n
    code += evaluatebasis.compute_scaling(element, Indent, format)

    # Compute auxilliary functions currently only 2D and 3D is supported
    if (element.cell_shape() == 2):
        code += evaluatebasis.compute_psitilde_a(element, Indent, format)
        code += evaluatebasis.compute_psitilde_b(element, Indent, format)
    elif (element.cell_shape() == 3):
        code += evaluatebasis.compute_psitilde_a(element, Indent, format)
        code += evaluatebasis.compute_psitilde_b(element, Indent, format)
        code += evaluatebasis.compute_psitilde_c(element, Indent, format)
    else:
        raise RuntimeError(), "Cannot compute auxilliary functions for shape: %d" % element.cell_shape()

    # Compute the basisvalues
    code += evaluatebasis.compute_basisvalues(element, Indent, format)

    return code

def tabulate_dmats(element, Indent, format):

    code = []

    # Get derivative matrices (coefficients) of basis functions, computed by FIAT at compile time
    derivative_matrices = element.basis().base.dmats

    # Get the shape of the element
    cell_shape = element.cell_shape()

    code += [Indent.indent(format["comment"]("Tables of derivative matrices (transpose)"))]

    # Generate tables for each spatial direction
    for i in range(cell_shape):

        # Extract derivatives for current direction (take transpose, FIAT ScalarPolynomialSet.deriv_all())
        matrix = numpy.transpose(derivative_matrices[i])

        # Get polynomial dimension of basis
        poly_dim = len(element.basis().base.bs)

        # Declare varable name for coefficients
        name = format["table declaration"] + "dmats%d[%d][%d]" %(i, poly_dim, poly_dim)
        value = evaluatebasis.tabulate_matrix(matrix, Indent, format)
        code += [(Indent.indent(name), Indent.indent(value))] + [""]

    return code

def compute_values(element, Indent, format):

    code = []

    # Get number of components, must change for tensor valued elements
    num_components = element.value_dimension(0)

    # Get polynomial dimension of basis
    poly_dim = len(element.basis().base.bs)

    # Get element shape
    cell_shape = element.cell_shape()

    code += [Indent.indent(format["comment"]("Compute value(s)"))]
    code += [Indent.indent(format["comment"]("Declare coefficients"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent(format["float declaration"] + "coeff%d_%d") % (i,j), "0.0")]
    code += [""]

    code += [Indent.indent(format["comment"]("Loop possible derivatives"))]
    code += [Indent.indent(format["loop"]("num_deriv", "n * %d" % cell_shape))]
    code += [Indent.indent(format["block begin"])]
    # Increase indentation
    Indent.increase()

    code += [Indent.indent(format["comment"]("Declare new coefficients"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent(format["float declaration"] + "new_coeff%d_%d") % (i,j),\
                                    "coefficients%d[dof][%d]" % (i,j))]

    # Debug coefficients
    value = "std::cout << "
    for i in range(num_components):
        for j in range(poly_dim):
            value += 'new_coeff%d_%d << " " << ' % (i,j)
    value += "std::endl;"
    code += [value]
    code += [""]

    # Update old coefficients
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent("coeff%d_%d") % (i,j), "new_coeff%d_%d" % (i,j))]
    code += [""]

#    code += [Indent.indent(format["loop"]("num_dir", "%d" % cell_shape))]
#    code += [Indent.indent(format["block begin"])]
#    # Increase indentation
#    Indent.increase()

    code += multiply_coeffs(element, Indent, format)
    code += dot_product(element, Indent, format)

    # Decrease indentation
    Indent.decrease()
    code += [Indent.indent(format["block end"])]

    return code

def multiply_coeffs(element, Indent, format):

    code = []

    # Get number of components, must change for tensor valued elements
    num_components = element.value_dimension(0)

    # Get polynomial dimension of basis
    poly_dim = len(element.basis().base.bs)

    # Get the shape of the element
    cell_shape = element.cell_shape()

    if (cell_shape == 2):

        for i in range(cell_shape):
            code += [Indent.indent("if(num_deriv == %d)\n{" % i)]
            # Increase indentation
            Indent.increase()

            for j in range(num_components):
                for k in range(poly_dim):
                    name = "new_coeff%d_%d" % (j,k)
                    value = format["add"]( [format["multiply"](["coeff%d_%d" %(j,l), "dmats%d[%d][%d]" %(i,l,k)]) \
                                            for l in range(poly_dim)])
                    code += [(Indent.indent(name), value)]


            # Debug coefficients
            value = "std::cout << "
            for i in range(num_components):
                for j in range(poly_dim):
                    value += 'new_coeff%d_%d << " " << ' % (i,j)
            value += "std::endl;"
            code += [value]

            # Decrease indentation
            Indent.decrease()
            code += [Indent.indent(format["block end"])]
    else:
        raise RuntimeError(), "Not implemented for 3D"

    return code + [""]

def dot_product(element, Indent, format):
    """This function computes the value of the basisfunction as the dot product of the
    coefficients and basisvalues """

    code = []

    # Get number of components, must change for tensor valued elements
    num_components = element.value_dimension(0)

    # Get polynomial dimension of basis
    poly_dim = len(element.basis().base.bs)

    for i in range(num_components):
        for j in range(poly_dim):
            name = "values[num_deriv]"
            value = format["add"]( [format["multiply"](["new_coeff%d_%d" %(i,k), "basisvalues[%d]" % k]) \
                                    for k in range(poly_dim)])
        code += [(Indent.indent(name), value)]

    return code



