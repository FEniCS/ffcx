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

# FFC code generation common modules
from utils import *

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

    # Compute number of derivatives that has to be computed, and declare an array to hold
    # the values of the derivatives on the reference element
    code += compute_num_derivatives(element, Indent, format)

    # Generate all possible combinations of derivatives
    code += generate_combinations(element, Indent, format)

    # Reset all values
    code += reset_values(element, Indent, format)

    # Check if we have just one element
    if (element.num_sub_elements() == 1):

        # Map degree of freedom to local degree of freedom for current element
        code += evaluatebasis.dof_map(0, Indent, format)

        code += generate_element_code(element, 0, False, Indent, format)

    # If the element is vector valued or mixed
    else:

        code += mixed_elements(element, Indent, format)

    return code

def compute_num_derivatives(element, Indent, format):
    """Computes the number of derivatives of order 'n' as: element.cell_shape()^n. Also
    declares a pointer to an array that holds the derivatives on the reference (FIAT) element."""

    code = []

    format_comment            = format["comment"]
    format_num_derivatives    = format["num derivatives"]
    format_float_declaration  = format["float declaration"]

    code += [format_comment("Compute number of derivatives")]

    code += [(Indent.indent(format["uint declaration"] + format_num_derivatives), "1")] + [""]

    # Loop order (n) to compute shape^n, std::pow doesn't work with (int, int) ambiguous call??
    code += [Indent.indent(format["loop"]("j", 0, format["argument derivative order"]))]

    # Increase indentation
    Indent.increase()

    code += [Indent.indent(format["times equal"](format_num_derivatives, element.cell_shape()))] + [""]

    # Decrease indentation
    Indent.decrease()

    # Declare pointer to array that holds derivatives on the FIAT element
    code += [format_comment("Declare pointer to array of derivatives on FIAT element")]

    code += [(Indent.indent(format_float_declaration + format["pointer"] + format["reference derivatives"]),\
              format["new"] + format_float_declaration + format["array access"](format_num_derivatives))]

    # Debug code
    code += [format["comment"]("Debug code")]
    code += [Indent.indent('std::cout << "number of derivatives = " << num_derivatives << std::endl;')]

    return code + [""]

def generate_combinations(element, Indent, format):
    "Generate all possible combinations of derivatives of order 'n'"

    code = []

    shape = element.cell_shape() - 1

    # Use code from codesnippets.py
    code += [Indent.indent(format["snippet combinations"])\
            % {"combinations": format["derivative combinations"], "shape-1": shape,\
               "num_derivatives" : format["num derivatives"], "n": format["argument derivative order"]}]

    # Debug code
    code += debug_combinations(element, Indent, format)
    
    return code + [""]

def generate_transform(element, Indent, format):
    """Generate the transformation matrix, whic is used to transform derivatives from reference
    element back to the physical element."""

    code = []

    # Generate code to construct the inverse of the Jacobian
    transform = "transform"
    num_derivatives = "num_derivatives"
    n = "n"
    combinations = "combinations"
    Jinv = "Jinv"

    if (element.cell_shape() == 2):

        # Use code from codesnippets.py
        code += [Indent.indent(format["snippet transform2D"])\
        % {"transform": format["transform matrix"], "num_derivatives" : format["num derivatives"],\
           "n": format["argument derivative order"], "combinations": format["derivative combinations"],\
           "Jinv":format["transform Jinv"]}]

    elif (element.cell_shape() == 3):
        code += [Indent.indent(format["snippet transform3D"])\
        % {"transform": format["transform matrix"], "num_derivatives" : format["num derivatives"],\
           "n": format["argument derivative order"], "combinations": format["derivative combinations"],\
           "Jinv":format["transform Jinv"]}]
    else:
        raise RuntimeError, "Cannot generate transform for shape: %d" %(element.cell_shape())

    # Debug code
    code += debug_transform(element, Indent, format)

    return code + [""]

def reset_values(element, Indent, format):
    "Reset all components of the 'values' array as it is a pointer to an array."

    code = []

    code += [Indent.indent(format["comment"]("Reset values"))]

    # Get number of components, change for tensor valued elements
    num_components = element.value_dimension(0)

    # Loop all values and set them equal to zero
    num_values = format["multiply"](["%d" %num_components, format["num derivatives"]])
    code += [Indent.indent(format["loop"]("j", 0, num_values))]

    # Increase indentation
    Indent.increase()

    # Reset values as it is a pointer
    code += [(Indent.indent(format["argument values"] + format["array access"]("j")),format["floating point"](0.0))]

    # Decrease indentation
    Indent.decrease()

    return code + [""]

def generate_element_code(element, value_num, vector, Indent, format):
    "Generate code for each basis element"

    code = []

    # Compute basisvalues, from evaluatebasis.py
    code += evaluatebasis.generate_basisvalues(element, Indent, format)

    # Tabulate coefficients
    code += evaluatebasis.tabulate_coefficients(element, Indent, format)

    code += [Indent.indent(format["comment"]("Interesting (new) part"))]

    # Tabulate coefficients for derivatives
    code += tabulate_dmats(element, Indent, format)

    # Compute the number of derivatives
    code += compute_num_derivatives(element, Indent, format)

    # Generate the combinations of derivatives
    code += generate_combinations(element, Indent, format)

    # Generate the transformation matrix
    code += generate_transform(element, Indent, format)

    # Compute the derivatives of the basisfunctions on the reference (FIAT) element, 
    # as the dot product of the new coefficients and basisvalues
    code += compute_reference_derivatives(element, Indent, format)

    # Transform derivatives to physical element by multiplication with the transformation matrix
    code += transform_derivatives(element, Indent, format)

    # Delete pointers
    code += delete_pointers(element, Indent, format)

    return code

def mixed_elements(element, Indent, format):
    "Generate code for each sub-element in the event of vector valued elements or mixed elements"

    code = []

    # Prefetch formats to speed up code generation
    format_block_begin = format["block begin"]
    format_block_end = format["block end"]
    format_dof_map_if = format["dof map if"]

    # Extract basis elements, and determine number of elements
    elements = evaluatebasis.extract_elements(element)
    num_elements = len(elements)

    sum_value_num = 0
    sum_space_dim = 0

    # Generate code for each element
    for i in range(num_elements):

        # Get sub element
        basis_element = elements[i]

        # FIXME: This must most likely change for tensor valued elements
        value_dimension = basis_element.value_dimension(0)
        value_num = basis_element.value_dimension(0)
        space_dim = basis_element.space_dimension()

        # Determine if the element has a value, for the given dof
        code += [Indent.indent(format_dof_map_if(sum_space_dim, sum_space_dim + space_dim -1))]
        code += [Indent.indent(format_block_begin)]
        # Increase indentation
        Indent.increase()

        # Generate map from global to local dof
        code += evaluatebasis.dof_map(sum_space_dim, Indent, format)

        # Generate code for basis element
        code += generate_element_code(basis_element, sum_value_num, True, Indent, format)

        # Decrease indentation, finish block - end element code
        Indent.decrease()
        code += [Indent.indent(format_block_end)] + [""]

        # Increase sum of value dimension, and space dimension
        sum_value_num += value_num
        sum_space_dim += space_dim

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

def compute_reference_derivatives(element, Indent, format):

    code = []

    # Get number of components, must change for tensor valued elements
    num_components = element.value_dimension(0)

    # Get polynomial dimension of basis
    poly_dim = len(element.basis().base.bs)

    # Get element shape
    cell_shape = element.cell_shape()

    num_derivatives = "num_derivatives"
    n = "n" # order of derivative
    derivatives = "derivatives"

    code += [Indent.indent(format["comment"]("Compute value(s)"))]
    code += [Indent.indent(format["comment"]("Declare coefficients"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent(format["float declaration"] + "coeff%d_%d") % (i,j), "0.0")]
    code += [""]

    code += [Indent.indent(format["comment"]("Declare new coefficients"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent(format["float declaration"] + "new_coeff%d_%d") % (i,j),\
                                    "coefficients%d[dof][%d]" % (i,j))]
    code += [""]

    code += [Indent.indent(format["comment"]("Loop possible derivatives"))]
    code += [Indent.indent(format["loop"]("num_deriv", 0, num_derivatives))]
    code += [Indent.indent(format["block begin"])]
    # Increase indentation
    Indent.increase()

    code += [Indent.indent(format["comment"]("Get values from coefficients array"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent("new_coeff%d_%d") % (i,j),\
                                    "coefficients%d[dof][%d]" % (i,j))]
    code += [""]


    # Debug coefficients
#    value = "std::cout << "
#    for i in range(num_components):
#        for j in range(poly_dim):
#            value += 'new_coeff%d_%d << " " << ' % (i,j)
#    value += "std::endl;"
#    code += [value]
#    code += [""]

    code += [Indent.indent(format["comment"]("Loop order"))]
    code += [Indent.indent(format["loop"]("j", 0, n))]
    code += [Indent.indent(format["block begin"])]
    # Increase indentation
    Indent.increase()

    # Update old coefficients
    code += [Indent.indent(format["comment"]("Update old coefficients"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent("coeff%d_%d") % (i,j), "new_coeff%d_%d" % (i,j))]
    code += [""]

    # Update new coefficients
    code += multiply_coeffs(element, Indent, format)

    # Decrease indentation
    Indent.decrease()
    code += [Indent.indent(format["block end"])]

    # Compute derivatives on reference element
    code += [Indent.indent(format["comment"]\
    ("Compute derivatives on reference element as dot product of coefficients and basisvalues"))]
    for i in range(num_components):
        for j in range(poly_dim):
            name = "%s[num_deriv]" %derivatives
            value = format["add"]( [format["multiply"](["new_coeff%d_%d" %(i,k), "basisvalues[%d]" % k]) \
                                    for k in range(poly_dim)])
        code += [(Indent.indent(name), value)]

    # Decrease indentation
    Indent.decrease()
    code += [Indent.indent(format["block end"])]

    # Debug code
    code += [Indent.indent('std::cout << "%s = " << std::endl;' %derivatives)]
    code += [Indent.indent(format["loop"]("j", 0, num_derivatives))]
    # Increase indent
    Indent.increase()
    code += [Indent.indent('std::cout << %s << " ";') % (derivatives + "[j]")]

    # Decrease indent
    Indent.decrease()
    code += [Indent.indent("std::cout << std::endl;")]

    return code + [""]

def multiply_coeffs(element, Indent, format):

    code = []

    # Get number of components, must change for tensor valued elements
    num_components = element.value_dimension(0)

    # Get polynomial dimension of basis
    poly_dim = len(element.basis().base.bs)

    # Get the shape of the element
    cell_shape = element.cell_shape()

    num_deriv = "num_deriv"
    combinations = "combinations"
    num_n = "j" # index for looping 'n'

    if (cell_shape == 2):

        for i in range(cell_shape):

# not language-independent
            code += [Indent.indent("if(%s[%s][%s] == %d)" % (combinations, num_deriv, num_n, i))]
            code += [Indent.indent(format["block begin"])]
            # Increase indentation
            Indent.increase()
            for j in range(num_components):
                for k in range(poly_dim):
                    name = "new_coeff%d_%d" % (j,k)
                    value = format["add"]( [format["multiply"](["coeff%d_%d" %(j,l), "dmats%d[%d][%d]" %(i,l,k)]) \
                                            for l in range(poly_dim)])
                    code += [(Indent.indent(name), value)]


            # Debug coefficients
#            value = "std::cout << "
#            for i in range(num_components):
#                for j in range(poly_dim):
#                    value += 'new_coeff%d_%d << " " << ' % (i,j)
#            value += "std::endl;"
#            code += [value]

            # Decrease indentation
            Indent.decrease()
            code += [Indent.indent(format["block end"])]
    else:
        raise RuntimeError, "Not implemented for 3D"

    return code + [""]

def transform_derivatives(element, Indent, format):
    """This function computes the value of the basisfunction as the dot product of the
    coefficients and basisvalues """

    code = []

    num_derivatives = "num_derivatives"
    values = "values"
    transform = "transform"
    derivatives = "derivatives"

    code += [Indent.indent(format["loop"]("j", 0, num_derivatives))]
    # Increase indentation
    Indent.increase()
    code += [(Indent.indent(values + "[j]"), "0.0")]
    # Decrease indentation
    Indent.decrease()
    
    code += [Indent.indent(format["loop"]("row", 0, num_derivatives))]
    code += [Indent.indent(format["block begin"])]
    # Increase indentation
    Indent.increase()
    code += [Indent.indent(format["loop"]("col", 0, num_derivatives))]
    # Increase indentation
    Indent.increase()

    # Compute values
    code += [Indent.indent(format["add equal"](values + "[row]",\
    format["multiply"]([transform + "[row][col]", derivatives + "[col]"])))]

    # Decrease indentation
    Indent.decrease()

    # Decrease indentation
    Indent.decrease()

    code += [Indent.indent(format["block end"])]

    return code



def delete_pointers(element, Indent, format):

    code = []

    # Delete pointers
    code += [format["comment"]("Delete pointer to array of derivatives on FIAT element")]
    code += [Indent.indent("delete [] derivatives;")] + [""]

    code += [format["comment"]("Delete pointer to array of combinations of derivatives")]
    code += [Indent.indent("delete [] combinations;")] + [""]

    return code + [""]



def debug_combinations(element, Indent, format):

    code = []

    # Debug code
    code += [Indent.indent('std::cout << "%s = " << std::endl;' % format["derivative combinations"])]
    code += [Indent.indent(format["loop"]("j", 0, format["num derivatives"]))]
    code += [Indent.indent(format["block begin"])]
    # Increase indent
    Indent.increase()
    code += [Indent.indent(format["loop"]("k", 0, format["argument derivative order"]))]
    code += [Indent.indent(format["block begin"])]
    # Increase indent
    Indent.increase()
    code += [Indent.indent('std::cout << %s << " ";') % (format["derivative combinations"] + "[j][k]")]

    # Decrease indent
    Indent.decrease()
    code += [Indent.indent(format["block end"])]
    code += [Indent.indent("std::cout << std::endl;")]

    # Decrease indent
    Indent.decrease()
    code += [Indent.indent(format["block end"])]

    return code

def debug_transform(element, Indent, format):

    code = []

    cell_shape = element.cell_shape()
    num_derivatives = format["num derivatives"]
    Jinv = format["transform Jinv"]
    transform = format["transform matrix"]

    # Debug code
    code += [format["comment"]("Debug code")]
    code += [Indent.indent("std::cout.precision(3);")]

    # Jinv
    code += [Indent.indent('std::cout << "%s = " << std::endl;' % Jinv)]
    code += [Indent.indent(format["loop"]("j", 0, cell_shape))]
    code += [Indent.indent(format["block begin"])]
    # Increase indent
    Indent.increase()
    code += [Indent.indent(format["loop"]("k", 0, cell_shape))]
    code += [Indent.indent(format["block begin"])]
    # Increase indent
    Indent.increase()
    code += [Indent.indent('std::cout << %s << " ";') % (Jinv + "[j][k]")]

    # Decrease indent
    Indent.decrease()
    code += [Indent.indent(format["block end"])]
    code += [Indent.indent("std::cout << std::endl;")]

    # Decrease indent
    Indent.decrease()
    code += [Indent.indent(format["block end"])]

    # Transform matrix
    code += [Indent.indent('std::cout << "%s = " << std::endl;' % transform)]
    code += [Indent.indent(format["loop"]("j", 0, num_derivatives))]
    code += [Indent.indent(format["block begin"])]
    # Increase indent
    Indent.increase()
    code += [Indent.indent(format["loop"]("k", 0, num_derivatives))]
    code += [Indent.indent(format["block begin"])]
    # Increase indent
    Indent.increase()
    code += [Indent.indent('std::cout << %s << " ";') % (transform + "[j][k]")]

    # Decrease indent
    Indent.decrease()
    code += [Indent.indent(format["block end"])]
    code += [Indent.indent("std::cout << std::endl;")]

    # Decrease indent
    Indent.decrease()
    code += [Indent.indent(format["block end"])]

    return code


