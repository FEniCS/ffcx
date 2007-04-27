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
# clean code

    code = []

    Indent = evaluatebasis.IndentControl()

    # Get coordinates and generate map
    code += evaluatebasis.generate_map(element, Indent, format)

    # Compute number of derivatives that has to be computed, and declare an array to hold
    # the values of the derivatives on the reference element
    code += compute_num_derivatives(element, Indent, format)

    # Generate all possible combinations of derivatives
    code += generate_combinations(element, Indent, format)

    # Generate the transformation matrix
    code += generate_transform(element, Indent, format)

    # Reset all values
    code += reset_values(element, Indent, format)

    # Check if we have just one element
    if (element.num_sub_elements() == 1):

        # Map degree of freedom to local degree of freedom for current element
        code += evaluatebasis.dof_map(0, Indent, format)

        code += generate_element_code(element, 0, Indent, format)

    # If the element is vector valued or mixed
    else:

        code += mixed_elements(element, Indent, format)

    return code

def compute_num_derivatives(element, Indent, format):
    "Computes the number of derivatives of order 'n' as: element.cell_shape()^n."

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

    # Debug code
#    code += [Indent.indent(format["comment"]("Debug code"))]
#    code += [Indent.indent('std::cout << "number of derivatives = " << num_derivatives << std::endl;')]

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
#    code += debug_combinations(element, Indent, format)
    
    return code + [""]

def generate_transform(element, Indent, format):
    """Generate the transformation matrix, whic is used to transform derivatives from reference
    element back to the physical element."""

    code = []

    # Generate code to construct the inverse of the Jacobian, use code from codesnippets.py
    # 2D
    if (element.cell_shape() == 2):
        code += [Indent.indent(format["snippet transform2D"])\
        % {"transform": format["transform matrix"], "num_derivatives" : format["num derivatives"],\
           "n": format["argument derivative order"], "combinations": format["derivative combinations"],\
           "Jinv":format["transform Jinv"]}]
    # 3D
    elif (element.cell_shape() == 3):
        code += [Indent.indent(format["snippet transform3D"])\
        % {"transform": format["transform matrix"], "num_derivatives" : format["num derivatives"],\
           "n": format["argument derivative order"], "combinations": format["derivative combinations"],\
           "Jinv":format["transform Jinv"]}]
    else:
        raise RuntimeError, "Cannot generate transform for shape: %d" %(element.cell_shape())

    # Debug code
#    code += debug_transform(element, Indent, format)

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

def generate_element_code(element, sum_value_dim, Indent, format):
    "Generate code for each basis element"

    code = []

    # Compute basisvalues, from evaluatebasis.py
    code += evaluatebasis.generate_basisvalues(element, Indent, format)

    # Tabulate coefficients
    code += evaluatebasis.tabulate_coefficients(element, Indent, format)

    code += [Indent.indent(format["comment"]("Interesting (new) part"))]

    # Tabulate coefficients for derivatives
    code += tabulate_dmats(element, Indent, format)

    # Compute the derivatives of the basisfunctions on the reference (FIAT) element, 
    # as the dot product of the new coefficients and basisvalues
    code += compute_reference_derivatives(element, Indent, format)

    # Transform derivatives to physical element by multiplication with the transformation matrix
    code += transform_derivatives(element, sum_value_dim, Indent, format)

    # Delete pointers
    code += delete_pointers(element, Indent, format)

    return code

def mixed_elements(element, Indent, format):
    "Generate code for each sub-element in the event of vector valued elements or mixed elements"

    code = []

    # Prefetch formats to speed up code generation
    format_block_begin  = format["block begin"]
    format_block_end    = format["block end"]
    format_dof_map_if   = format["dof map if"]

    # Extract basis elements, and determine number of elements
    elements = evaluatebasis.extract_elements(element)
    num_elements = len(elements)

    sum_value_dim = 0
    sum_space_dim = 0

    # Generate code for each element
    for i in range(num_elements):

        # Get sub element
        basis_element = elements[i]

        # FIXME: This must most likely change for tensor valued elements
        value_dim = basis_element.value_dimension(0)
        space_dim = basis_element.space_dimension()

        # Determine if the element has a value, for the given dof
        code += [Indent.indent(format_dof_map_if(sum_space_dim, sum_space_dim + space_dim -1))]
        code += [Indent.indent(format_block_begin)]
        # Increase indentation
        Indent.increase()

        # Generate map from global to local dof
        code += evaluatebasis.dof_map(sum_space_dim, Indent, format)

        # Generate code for basis element
        code += generate_element_code(basis_element, sum_value_dim, Indent, format)

        # Decrease indentation, finish block - end element code
        Indent.decrease()
        code += [Indent.indent(format_block_end)] + [""]

        # Increase sum of value dimension, and space dimension
        sum_value_dim += value_dim
        sum_space_dim += space_dim

    return code

def tabulate_dmats(element, Indent, format):
    "Tabulate the derivatives of the polynomial base"

    code = []

    # Prefetch formats to speed up code generation
    format_table          = format["table declaration"]
    format_dmats          = format["dmats table"]
    format_matrix_access  = format["matrix access"]

    # Get derivative matrices (coefficients) of basis functions, computed by FIAT at compile time
    derivative_matrices = element.basis().base.dmats

    # Get the shape of the element
    cell_shape = element.cell_shape()

    code += [Indent.indent(format["comment"]("Tables of derivatives of the polynomial base (transpose)"))]

    # Generate tables for each spatial direction
    for i in range(cell_shape):

        # Extract derivatives for current direction (take transpose, FIAT ScalarPolynomialSet.deriv_all())
        matrix = numpy.transpose(derivative_matrices[i])

        # Get polynomial dimension of basis
        poly_dim = len(element.basis().base.bs)

        # Declare varable name for coefficients
        name = format_table + format_dmats(i) + format_matrix_access(poly_dim, poly_dim)
        value = tabulate_matrix(matrix, format)
        code += [(Indent.indent(name), Indent.indent(value))] + [""]

    return code

def compute_reference_derivatives(element, Indent, format):
    """Compute derivatives on the reference element by recursively multiply coefficients with
    the relevant derivatives of the polynomial base until the requested order of derivatives
    has been reached. After this take the dot product with the basisvalues."""

    code = []

    # Prefetch formats to speed up code generation
    format_comment          = format["comment"]
    format_float            = format["float declaration"]
    format_coeff            = format["coefficient scalar"]
    format_new_coeff        = format["new coefficient scalar"]
    format_secondary_index  = format["secondary index"]
    format_floating_point   = format["floating point"]
    format_num_derivatives  = format["num derivatives"]
    format_loop             = format["loop"]
    format_block_begin      = format["block begin"]
    format_block_end        = format["block end"]
    format_coefficients     = format["coefficients table"]
    format_dof              = format["local dof"]
    format_n                = format["argument derivative order"]
    format_derivatives      = format["reference derivatives"]
    format_matrix_access    = format["matrix access"]
    format_array_access     = format["array access"]
    format_add              = format["add"]
    format_multiply         = format["multiply"]
    format_basisvalue       = format["basisvalue"]

    # Get number of components, must change for tensor valued elements
    num_components = element.value_dimension(0)

    # Get polynomial dimension of basis
    poly_dim = len(element.basis().base.bs)

    # Get element shape
    cell_shape = element.cell_shape()

    code += [Indent.indent(format_comment("Compute reference derivatives"))]

    # Declare pointer to array that holds derivatives on the FIAT element
    code += [Indent.indent(format_comment("Declare pointer to array of derivatives on FIAT element"))]

    # The size of the array of reference derivatives is equal to the number of derivatives
    # times the value dimension of the basis element
    if (num_components == 1):
        code += [(Indent.indent(format_float + format["pointer"] + format["reference derivatives"]),\
                  format["new"] + format_float + format["array access"](format_num_derivatives))]
    else:
        code += [(Indent.indent(format_float + format["pointer"] + format["reference derivatives"]),\
                  format["new"] + format_float + format["array access"]\
                  (format_multiply(["%s" %num_components, format_num_derivatives])) )]

    code += [""]

    code += [Indent.indent(format_comment("Declare coefficients"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent(format_float + format_coeff(i) + format_secondary_index(j)),\
                      format_floating_point(0.0))]
    code += [""]

    code += [Indent.indent(format_comment("Declare new coefficients"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent(format_float + format_new_coeff(i) + format_secondary_index(j)),\
                      format_floating_point(0.0))]
    code += [""]

    code += [Indent.indent(format_comment("Loop possible derivatives"))]
    code += [Indent.indent(format_loop("deriv_num", 0, format_num_derivatives))]
    code += [Indent.indent(format_block_begin)]
    # Increase indentation
    Indent.increase()

    code += [Indent.indent(format_comment("Get values from coefficients array"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent(format_new_coeff(i) + format_secondary_index(j)),\
                      format_coefficients(i) + format_matrix_access(format_dof, j))]
    code += [""]

    code += [Indent.indent(format_comment("Loop derivative order"))]
    code += [Indent.indent(format_loop("j", 0, format_n))]
    code += [Indent.indent(format_block_begin)]
    # Increase indentation
    Indent.increase()

    # Update old coefficients
    code += [Indent.indent(format_comment("Update old coefficients"))]
    for i in range(num_components):
        for j in range(poly_dim):
            code += [(Indent.indent(format_coeff(i) + format_secondary_index(j)),\
                      format_new_coeff(i) + format_secondary_index(j))]
    code += [""]

    # Update new coefficients
    code += multiply_coeffs(element, Indent, format)

    # Decrease indentation
    Indent.decrease()
    code += [Indent.indent(format_block_end)]

    # Compute derivatives on reference element
    code += [Indent.indent(format_comment\
    ("Compute derivatives on reference element as dot product of coefficients and basisvalues"))]
    for i in range(num_components):
        if (i == 0):
            name = format_derivatives + format_array_access("deriv_num")
        elif (i == 1):
            name = format_derivatives + format_array_access(format_add([format_num_derivatives, "deriv_num"] ))
        else:
            name = format_derivatives + format_array_access(format_add([format_multiply\
                   (["%d" %i, format_num_derivatives]), "deriv_num"] ))

        value = format_add([ format_multiply([format_new_coeff(i) + format_secondary_index(k),\
                             format_basisvalue(k)]) for k in range(poly_dim) ])
        code += [(Indent.indent(name), value)]

    # Decrease indentation
    Indent.decrease()
    code += [Indent.indent(format_block_end)]

    return code + [""]

def multiply_coeffs(element, Indent, format):
    "Auxilliary function that multiplies coefficients with directional derivatives."

    code = []

    # Prefetch formats to speed up code generation
    format_if               = format["if"]
    format_group            = format["grouping"]
    format_isequal          = format["is equal"]
    format_combinations     = format["derivative combinations"]
    format_block_begin      = format["block begin"]
    format_block_end        = format["block end"]
    format_coeff            = format["coefficient scalar"]
    format_new_coeff        = format["new coefficient scalar"]
    format_secondary_index  = format["secondary index"]
    format_add              = format["add"]
    format_multiply         = format["multiply"]
    format_dmats            = format["dmats table"]
    format_matrix_access    = format["matrix access"]

    # Get number of components, must change for tensor valued elements
    num_components = element.value_dimension(0)

    # Get polynomial dimension of basis
    poly_dim = len(element.basis().base.bs)

    # Get the shape of the element
    cell_shape = element.cell_shape()

    for i in range(cell_shape):

# not language-independent
        code += [Indent.indent(format_if + format_group( format_combinations + \
                 format_matrix_access("deriv_num","j") + format_isequal + "%d" %(i) ) )]
        code += [Indent.indent(format_block_begin)]
        # Increase indentation
        Indent.increase()
        for j in range(num_components):
            for k in range(poly_dim):
                name = format_new_coeff(j) + format_secondary_index(k)
                value = format_add( [format_multiply([format_coeff(j) + format_secondary_index(l),\
                        format_dmats(i) + format_matrix_access(l,k)]) for l in range(poly_dim)])
                code += [(Indent.indent(name), value)]

        # Decrease indentation
        Indent.decrease()
        code += [Indent.indent(format_block_end)]

    return code + [""]

def transform_derivatives(element, sum_value_dim, Indent, format):
    """This function computes the value of the basisfunction as the dot product of the
    coefficients and basisvalues """

    code = []

    # Prefetch formats to speed up code generation
    format_loop             = format["loop"]
    format_num_derivatives  = format["num derivatives"]
    format_derivatives      = format["reference derivatives"]
    format_values           = format["argument values"] 
    format_multiply         = format["multiply"]
    format_add              = format["add"]
    format_array_access     = format["array access"]
    format_matrix_access    = format["matrix access"]
    format_transform        = format["transform matrix"]

    code += [Indent.indent(format["comment"]("Transform derivatives back to physical element"))]
    
    code += [Indent.indent(format_loop("row", 0, format_num_derivatives))]
    code += [Indent.indent(format["block begin"])]
    # Increase indentation
    Indent.increase()
    code += [Indent.indent(format_loop("col", 0, format_num_derivatives))]
    code += [Indent.indent(format["block begin"])]
    # Increase indentation
    Indent.increase()

    # Get number of components, must change for tensor valued elements
    num_components = element.value_dimension(0)

    # Compute offset in array values if any
    for i in range(num_components):
        if (sum_value_dim + i == 0):
            name = format_values + format_array_access("row")
        elif (sum_value_dim + i == 1):
            name = format_values + format_array_access(format_add([format_num_derivatives, "row"]))
        else:
            offset_name = format_multiply(["%d" %(sum_value_dim + i), format_num_derivatives])
            name = format_values + format_array_access(format_add([offset_name, "row"]))
        if (i == 0):
            value = format_multiply([format_transform + format_matrix_access("row","col"),\
                    format_derivatives + format_array_access("col")])
        elif (i == 1):
            value = format_multiply([format_transform + format_matrix_access("row","col"),\
                    format_derivatives + format_array_access(format_add([format_num_derivatives, "col"]))])
        else:
            offset_value = format_multiply(["%d" %i, format_num_derivatives])

            value = format_multiply([format_transform + format_matrix_access("row","col"),\
                    format_derivatives + format_array_access(format_add([offset_value, "col"]))])

        # Compute values
        code += [Indent.indent(format["add equal"](name, value))]

    # Decrease indentation
    Indent.decrease()
    code += [Indent.indent(format["block end"])]

    # Decrease indentation
    Indent.decrease()

    code += [Indent.indent(format["block end"])]

    return code

def delete_pointers(element, Indent, format):
    "Delete the pointers to arrays."

    code = []

    # Delete pointers
    code += [Indent.indent(format["comment"]("Delete pointer to array of derivatives on FIAT element"))]
    code += [Indent.indent(format["delete"] + format["array access"]("") + " " +\
                           format["reference derivatives"] + format["end line"])] + [""]

    code += [Indent.indent(format["comment"]("Delete pointer to array of combinations of derivatives"))]
    code += [Indent.indent(format["delete"] + format["array access"]("") + " " +\
                           format["derivative combinations"] + format["end line"])] + [""]

    return code


def debug_combinations(element, Indent, format):

    code = []

    code += [format["comment"]("Debug code")]
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

def debug_reference_derivatives(element, Indent, format):

    # Debug code
    code = [Indent.indent('std::cout << "%s = " << std::endl;' %format["reference derivatives"])]
    code += [Indent.indent(format["loop"]("j", 0, format["num derivatives"]))]
    # Increase indent
    Indent.increase()
    code += [Indent.indent('std::cout << %s << " ";') % (format["reference derivatives"] + "[j]")]

    # Decrease indent
    Indent.decrease()
    code += [Indent.indent("std::cout << std::endl;")]
    return code

def debug_():
    "misc"

    code = []

    # Debug coefficients
#    value = "std::cout << "
#    for i in range(num_components):
#        for j in range(poly_dim):
#            value += 'new_coeff%d_%d << " " << ' % (i,j)
#    value += "std::endl;"
#    code += [value]
#    code += [""]

            # Debug coefficients
#            value = "std::cout << "
#            for i in range(num_components):
#                for j in range(poly_dim):
#                    value += 'new_coeff%d_%d << " " << ' % (i,j)
#            value += "std::endl;"
#            code += [value]

    return code
