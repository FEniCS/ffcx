"""Code generation for evaluation of finite element basis values. This module generates
code which is more or less a C++ representation of FIAT code."""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-12-14"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-22

# Python modules
import math
import numpy

# FFC modules
from log import error, debug_code
from cpp import remove_unused
from codegeneratorsutils import IndentControl
from codegeneratorsutils import inner_product, tabulate_matrix, tabulate_vector
from quadrature.quadraturegenerator_utils import generate_loop

# Temporary import
from cpp import format_old as format

def _evaluate_basis(data):
    """Generate run time code to evaluate an element basisfunction at an
    arbitrary point. The value(s) of the basisfunction is/are
    computed as in FIAT as the dot product of the coefficients (computed at compile time)
    and basisvalues which are dependent on the coordinate and thus have to be computed at
    run time.

    Currently the following elements are supported in 1D:

    Lagrange                + mixed/vector valued
    Discontinuous Lagrange  + mixed/vector valued

    Currently the following elements are supported in 2D and 3D:

    Lagrange                + mixed/vector valued
    Discontinuous Lagrange  + mixed/vector valued
    Crouzeix-Raviart        + mixed/vector valued
    Brezzi-Douglas-Marini   + mixed/vector valued

    Not supported in 2D or 3D:

    Raviart-Thomas ? (not tested since it is broken in FFC, but should work)
    Nedelec (broken?)

    Tensor valued elements!"""
    if data is None:
        return ""
    code = []

    Indent = IndentControl()

    # Get coordinates and generate map
    # FIXME: KBO: If this remains simple, inline function here.
    code += _generate_map(data, Indent, format)

    # Check if we have just one element

    if (data["num_sub_elements"] == 1):

        # Reset values, change for tensor valued elements
        # FIXME: KBO: If this remains simple, inline function here.
        value_shape = data["value_shape"]
        code += _reset_values(value_shape, False, Indent, format)

        # Map degree of freedom to local degree of freedom for current element
        code += _map_dof(0, Indent, format)

        # Generate element code
        code += _generate_element_code(data, 0, False, Indent, format)

    # If the element is of type VectorElement or MixedElement
    else:
        # FIXME: KBO: Only support scalar elements for now
        error("Only scalar elements are supported.")
#        # Reset values, change for tensor valued elements
        # FIXME: KBO: If this remains simple, inline function here.
#        code += _reset_values(data.value_dimension(0), True, Indent, format)

#        # Generate element code, for all sub-elements
#        code += _mixed_elements(data, Indent, format)

    lines = format["generate body"](code)

    # Remove unused variables (from transformations and mappings) in code.
    code = remove_unused(lines)
    return code

def _generate_map(data, Indent, format):
    """Generates map from physical element to the UFC reference element, and from this element
    to reference square/cube. An affine map is assumed for the first mapping and for the second
    map this function implements the UFC version of the FIAT functions, eta_triangle( xi )
    and eta_tetrahedron( xi ) from reference.py"""

    code = []

    # Prefetch formats to speed up code generation
    format_comment        = format["comment"]
    format_floating_point = format["floating point"]
    format_epsilon        = format["epsilon"]

    cell_domain = "triangle"
#    # Get coordinates and map to the UFC reference element from codesnippets.py
#    code += [Indent.indent(format["coordinate map FIAT"](data.cell_domain()))] + [""]



    # Get coordinates and map to the FIAT reference element from codesnippets.py
    code += [Indent.indent(format["coordinate map FIAT"](cell_domain))] + [""]

    # FIXME: Verify that we don't need to apply additional transformations once we're on the reference element
#    if (data.cell_domain() == "interval"):

#        # Map coordinates to the reference interval
#        code += [Indent.indent(format_comment("Map coordinates to the reference interval"))]

#        # Code snippet reproduced from FIAT: reference.py: eta_line(xi)
#        code += [Indent.indent(format["snippet eta_interval"])]

#    elif (data.cell_domain() == "triangle"):

#        # Map coordinates to the reference square
#        code += [Indent.indent(format_comment("Map coordinates to the reference square"))]

#        # Code snippet reproduced from FIAT: reference.py: eta_triangle(xi)
#        code += [Indent.indent(format["snippet eta_triangle"]) %(format_floating_point(format_epsilon))]

#    elif (data.cell_domain() == "tetrahedron"):

#        # Map coordinates to the reference cube
#        code += [Indent.indent(format_comment("Map coordinates to the reference cube"))]

#        # Code snippet reproduced from FIAT: reference.py: eta_tetrahedron(xi)
#        code += [Indent.indent(format["snippet eta_tetrahedron"]) %(format_floating_point(format_epsilon),\
#                       format_floating_point(format_epsilon))]
#    else:
#        error("Cannot generate map for shape: %d" %(data.cell_domain()))

    return code

def _reset_values(num_components, vector, Indent, format):
    "Reset all components of the 'values' array as it is a pointer to an array."

    code = []

    if (vector or num_components != ()):
        # Reset values as it is a pointer
        code += [Indent.indent(format["comment"]("Reset values"))]
        code += [(Indent.indent(format["argument values"] + format["array access"](i)),\
                                format["floating point"](0.0)) for i in range(num_components)]
    else:
        code += [Indent.indent(format["comment"]("Reset values"))]
        code += [(Indent.indent(format["pointer"] + format["argument values"]), format["floating point"](0.0))]

    return code + [""]

def _map_dof(sum_space_dim, Indent, format):
    """This function creates code to map a basis function to a local basis function.
    Example, the following mixed element:

    element = VectorElement("Lagrange", "triangle", 2)

    has the element list, elements = [Lagrange order 2, Lagrange order 2] and 12 dofs (6 each).

    The evaluation of basis function 8 is then mapped to 2 (8-6) for local element no. 2."""

    # In case of only one element or the first element in a series then we don't subtract anything
    if sum_space_dim == 0:
        code = [Indent.indent(format["comment"]("Map degree of freedom to element degree of freedom"))]
        code += [(Indent.indent(format["const uint declaration"] + format["local dof"]), format["argument basis num"])]
    else:
        code = [Indent.indent(format["comment"]("Map degree of freedom to element degree of freedom"))]
        code += [(Indent.indent(format["const uint declaration"] + format["local dof"]),\
                format["subtract"]([format["argument basis num"], "%d" %sum_space_dim]))]

    return code + [""]

def _mixed_elements(data, Indent, format):
    "Generate code for each sub-element in the event of vector valued elements or mixed elements"

    code = []

    # Prefetch formats to speed up code generation
    format_block_begin = format["block begin"]
    format_block_end = format["block end"]
    format_dof_map_if = format["dof map if"]

    # Extract basis elements, and determine number of elements
    elements = data.extract_elements()
    num_elements = len(elements)

    sum_value_dim = 0
    sum_space_dim = 0

    # Generate code for each element
    for i in range(num_elements):

        # Get sub element
        basis_element = elements[i]

        # Get value and space dimension
        value_dim = basis_element.value_shape()[0]
        space_dim = basis_element.space_dimension()

        # Determine if the element has a value, for the given dof
        code += [Indent.indent(format_dof_map_if(sum_space_dim, sum_space_dim + space_dim -1))]
        code += [Indent.indent(format_block_begin)]
        # Increase indentation
        Indent.increase()

        # Generate map from global to local dof
        code += _map_dof(sum_space_dim, Indent, format)

        # Generate code for basis element
        code += _generate_element_code(basis_element, sum_value_dim, True, Indent, format)

        # Decrease indentation, finish block - end element code
        Indent.decrease()
        code += [Indent.indent(format_block_end)] + [""]

        # Increase sum of value dimension, and space dimension
        sum_value_dim += value_dim
        sum_space_dim += space_dim

    return code

def _generate_element_code(data, sum_value_dim, vector, Indent, format):
    "Generate code for a single basis element as the dot product of coefficients and basisvalues."

    code = []

    # Generate basisvalues
    code += _generate_basisvalues(data, Indent, format)

    # Tabulate coefficients
    code += _tabulate_coefficients(data, Indent, format)

    # Compute the value of the basisfunction as the dot product of the coefficients
    # and basisvalues
    code += _compute_values(data, sum_value_dim, vector, Indent, format)

    return code

def _generate_basisvalues(data, Indent, format):
    "Generate code to compute basis values"

    code = []

# FIXME: KBO: Scalings should not be needed anymore, delete
#    # Compute scaling of y and z 1/2(1-y)^n and 1/2(1-z)^n
#    code += _compute_scaling(data, Indent, format)

# FIXME: KBO: Like scalings this should not be needed anymore, delete
#    # Compute auxilliary functions
#    if data.cell_domain() == "interval":
#        code += _compute_psitilde_a(data, Indent, format)
#    elif data.cell_domain() == "triangle":
#        code += _compute_psitilde_a(data, Indent, format)
#        code += _compute_psitilde_b(data, Indent, format)
#    elif data.cell_domain() == "tetrahedron":
#        code += _compute_psitilde_a(data, Indent, format)
#        code += _compute_psitilde_b(data, Indent, format)
#        code += _compute_psitilde_c(data, Indent, format)
#    else:
#        error("Cannot compute auxilliary functions for shape: %d" %(data.cell_domain()))

# FIXME: KBO: If the above can be deleted we should remove this function and move
# code here.
    # Compute the basisvalues
    code += _compute_basisvalues(data, Indent, format)

    return code

def _tabulate_coefficients(data, Indent, format):
    """This function tabulates the element coefficients that are generated by FIAT at
    compile time."""

    code = []

    # Prefetch formats to speed up code generation
    format_comment            = format["comment"]
    format_table_declaration  = format["table declaration"]
    format_coefficients       = format["coefficients table"]
    format_matrix_access      = format["matrix access"]
    format_const_float        = format["const float declaration"]

    # Get coefficients from basis functions, computed by FIAT at compile time
    coefficients = data["coeffs"]

    # Scalar valued basis element [Lagrange, Discontinuous Lagrange, Crouzeix-Raviart]
    rank = data["value_rank"]
    if (rank == 0):
        coefficients = [coefficients]

    # Vector valued basis element [Raviart-Thomas, Brezzi-Douglas-Marini (BDM)]
    # FIXME: KBO: Verify that coefficients still look this way
    elif (rank == 1):
        error("Broken!")
        coefficients = numpy.transpose(coefficients, [1,0,2])
    else:
        error("Tensor elements not supported!")

    # Generate tables for each component
    code += [Indent.indent(format_comment("Table(s) of coefficients"))]
    for i, coeffs in enumerate(coefficients):

        # Get number of dofs and number of members of the expansion set
        num_dofs, num_mem = numpy.shape(coeffs)

        # Declare varable name for coefficients
        name = format_table_declaration + format_coefficients(i) + format_matrix_access(num_dofs, num_mem)
        value = tabulate_matrix(coeffs, format)

        # Generate array of values
        code += [(Indent.indent(name), Indent.indent(value))] + [""]

    return code

def _compute_values(data, sum_value_dim, vector, Indent, format):
    """This function computes the value of the basisfunction as the dot product of the
    coefficients and basisvalues """

    code = []

    # Prefetch formats to speed up code generation
    format_values           = format["argument values"]
    format_array_access     = format["array access"]
    format_matrix_access     = format["matrix access"]
#    format_add              = format["add"]

    print format

    format_multiply         = format["multiply"]
    format_coefficients     = format["coefficients table"]
    format_basisvalues      = format["basisvalues table"]
    format_j                = format["first free index"]
    format_dof              = format["local dof"]
#    format_secondary_index  = format["secondary index"]
#    format_basisvalue       = format["basisvalue"]
    format_pointer          = format["pointer"]
#    format_det              = format["determinant"]
#    format_inv              = format["inverse"]
#    format_add              = format["add"]
#    format_mult             = format["multiply"]
#    format_group            = format["grouping"]
#    format_tmp              = format["tmp declaration"]
#    format_tmp_access       = format["tmp access"]

    code += [Indent.indent(format["comment"]("Compute value(s)"))]

    # Get number of components, change for tensor valued elements
    shape = data["value_shape"]
    if len(shape) > 0:
        num_components = shape[0]
    else:
        num_components = 1

    # Check which transform we should use to map the basis functions
    # FIXME: KBO: Only support affine mapping for now
#    mapping = data.mapping()
#    if mapping == CONTRAVARIANT_PIOLA:
#        code += [Indent.indent(format["comment"]("Using contravariant Piola transform to map values back to the physical element"))]
#    elif mapping == COVARIANT_PIOLA:
#        code += [Indent.indent(format["comment"]("Using covariant Piola transform to map values back to the physical element"))]

    if (vector or num_components != 1):
        # Loop number of components
        for i in range(num_components):
            name = format_values + format_array_access(i + sum_value_dim)
#            value = format_add([format_multiply([format_coefficients(i) + format_secondary_index(j),\
#                    format_basisvalue(j)]) for j in range(poly_dim)])

#            # Use Piola transform to map basisfunctions back to physical data if needed
#            if mapping == CONTRAVARIANT_PIOLA:
#                code.insert(i+1,(Indent.indent(format_tmp(0, i)), value))
#                basis_col = [format_tmp_access(0, j) for j in range(data.cell().topological_dimension())]
#                jacobian_row = [format["transform"]("J", j, i, None) for j in range(data.cell().topological_dimension())]
#                inner = [format_mult([jacobian_row[j], basis_col[j]]) for j in range(data.cell().topological_dimension())]
#                sum = format_group(format_add(inner))
#                value = format_mult([format_inv(format_det(None)), sum])
#            elif mapping == COVARIANT_PIOLA:
#                code.insert(i+1,(Indent.indent(format_tmp(0, i)), value))
#                basis_col = [format_tmp_access(0, j) for j in range(data.cell().topological_dimension())]
#                inverse_jacobian_column = [format["transform"]("JINV", j, i, None) for j in range(data.cell().topological_dimension())]
#                inner = [format_mult([inverse_jacobian_column[j], basis_col[j]]) for j in range(data.cell().topological_dimension())]
#                sum = format_group(format_add(inner))
#                value = format_mult([sum])

#            code += [(Indent.indent(name), value)]
    else:
        # Get number of members of the expansion set
        # FIXME: KBO: Might be able to get this more elegantly, maybe move to
        # data or simply take the second value of the shape of the coefficients
        num_mem = data["num_expansion_members"]
        loop_vars = [(format_j, 0, num_mem)]

        name = format_pointer + format_values

        value = format_multiply([format_coefficients(0) + format_matrix_access(format_dof, format_j),\
                format_basisvalues + format_array_access(format_j)])

        line = [format["add equal"](name, value)]

        code += generate_loop(line, loop_vars, Indent, format)

    return code

# FIXME: KBO: Scalings should not be needed anymore, delete
#def _compute_scaling(data, Indent, format):
#    """Generate the scalings of y and z coordinates. This function is an implementation of
#    the FIAT function make_scalings( etas ) from expansions.py"""

#    code = []

#    # Prefetch formats to speed up code generation
#    format_y            = format["y coordinate"]
#    format_z            = format["z coordinate"]
#    format_grouping     = format["grouping"]
#    format_subtract     = format["subtract"]
#    format_multiply     = format["multiply"]
#    format_const_float  = format["const float declaration"]
#    format_scalings     = format["scalings"]

#    # Get the element degree
#    degree = data.degree()

#    # Get the element cell domain
#    element_cell_domain = data.cell_domain()

#    # For 1D scalings are not needed
#    if element_cell_domain == "interval":
#        code += [Indent.indent(format["comment"]("Generate scalings not needed for 1D"))]
#        return code + [""]
#    elif element_cell_domain == "triangle":
#        scalings = [format_y]
#        # Scale factor, for triangles 1/2*(1-y)^i i being the order of the element
#        factors = [format_grouping(format_subtract(["0.5", format_multiply(["0.5", format_y])]))]
#    elif element_cell_domain == "tetrahedron":
#        scalings = [format_y, format_z]
#        factors = [format_grouping(format_subtract(["0.5", format_multiply(["0.5", format_y])])),\
#                   format_grouping(format_subtract(["0.5", format_multiply(["0.5", format_z])]))]
#    else:
#        error("Cannot compute scaling for shape: %d" %(element_cell_domain))

#    code += [Indent.indent(format["comment"]("Generate scalings"))]

#    # Can be optimized by leaving out the 1.0 variable
#    for i in range(len(scalings)):

#      name = format_const_float + format_scalings(scalings[i], 0)
#      value = format["floating point"](1.0)
#      code += [(Indent.indent(name), value)]

#      for j in range(1, degree+1):
#          name = format_const_float + format_scalings(scalings[i], j)
#          value = format_multiply([format_scalings(scalings[i],j-1), factors[i]])
#          code += [(Indent.indent(name), value)]

#    return code + [""]

#def _compute_psitilde_a(data, Indent, format):
#    """Compute Legendre functions in x-direction. The function relies on
#    eval_jacobi_batch(a,b,n) to compute the coefficients.

#    The format is:
#    psitilde_a[0] = 1.0
#    psitilde_a[1] = a + b * x
#    psitilde_a[n] = a * psitilde_a[n-1] + b * psitilde_a[n-1] * x + c * psitilde_a[n-2]
#    where a, b and c are coefficients computed by eval_jacobi_batch(0,0,n)
#    and n is the element degree"""

#    code = []

#    # Get the element degree
#    degree = data.degree()

#    code += [Indent.indent(format["comment"]("Compute psitilde_a"))]

#    # Create list of variable names
#    variables = [format["x coordinate"], format["psitilde_a"]]

#    for n in range(degree+1):
#        # Declare variable
#        name = format["const float declaration"] + variables[1] + format["secondary index"](n)

#        # Compute value
#        value = _eval_jacobi_batch_scalar(0, 0, n, variables, format)

#        code += [(Indent.indent(name), value)]

#    return code + [""]

#def _compute_psitilde_b(data, Indent, format):
#    """Compute Legendre functions in y-direction. The function relies on
#    eval_jacobi_batch(a,b,n) to compute the coefficients.

#    The format is:
#    psitilde_bs_0[0] = 1.0
#    psitilde_bs_0[1] = a + b * y
#    psitilde_bs_0[n] = a * psitilde_bs_0[n-1] + b * psitilde_bs_0[n-1] * x + c * psitilde_bs_0[n-2]
#    psitilde_bs_(n-1)[0] = 1.0
#    psitilde_bs_(n-1)[1] = a + b * y
#    psitilde_bs_n[0] = 1.0
#    where a, b and c are coefficients computed by eval_jacobi_batch(2*i+1,0,n-i) with i in range(0,n+1)
#    and n is the element degree + 1"""

#    code = []

#    # Prefetch formats to speed up code generation
#    format_y                = format["y coordinate"]
#    format_psitilde_bs      = format["psitilde_bs"]
#    format_const_float      = format["const float declaration"]
#    format_secondary_index  = format["secondary index"]

#    # Get the element degree
#    degree = data.degree()

#    code += [Indent.indent(format["comment"]("Compute psitilde_bs"))]

#    for i in range(0, degree + 1):

#        # Compute constants for jacobi function
#        a = 2*i+1
#        b = 0
#        n = degree - i

#        # Create list of variable names
#        variables = [format_y, format_psitilde_bs(i)]

#        for j in range(n+1):
#            # Declare variable
#            name = format_const_float + variables[1] + format_secondary_index(j)

#            # Compute values
#            value = _eval_jacobi_batch_scalar(a, b, j, variables, format)

#            code += [(Indent.indent(name), value)]

#    return code + [""]

#def _compute_psitilde_c(data, Indent, format):
#    """Compute Legendre functions in y-direction. The function relies on
#    eval_jacobi_batch(a,b,n) to compute the coefficients.

#    The format is:
#    psitilde_cs_0[0] = 1.0
#    psitilde_cs_0[1] = a + b * y
#    psitilde_cs_0[n] = a * psitilde_cs_0[n-1] + b * psitilde_cs_0[n-1] * x + c * psitilde_cs_0[n-2]
#    psitilde_cs_(n-1)[0] = 1.0
#    psitilde_cs_(n-1)[1] = a + b * y
#    psitilde_cs_n[0] = 1.0
#    where a, b and c are coefficients computed by

#    [[jacobi.eval_jacobi_batch(2*(i+j+1),0, n-i-j) for j in range(0,n+1-i)] for i in range(0,n+1)]"""

#    code = []

#    # Prefetch formats to speed up code generation
#    format_z                = format["z coordinate"]
#    format_psitilde_cs      = format["psitilde_cs"]
#    format_const_float      = format["const float declaration"]
#    format_secondary_index  = format["secondary index"]

#    # Get the element degree
#    degree = data.degree()

#    code += [Indent.indent(format["comment"]("Compute psitilde_cs"))]

#    for i in range(0, degree + 1):
#        for j in range(0, degree + 1 - i):

#            # Compute constants for jacobi function
#            a = 2*(i+j+1)
#            b = 0
#            n = degree - i - j

#            # Create list of variable names
#            variables = [format_z, format_psitilde_cs(i,j)]

#            for k in range(n+1):
#                # Declare variable
#                name = format_const_float + variables[1] + format_secondary_index(k)

#                # Compute values
#                value = _eval_jacobi_batch_scalar(a, b, k, variables, format)

#                code += [(Indent.indent(name), value)]

#    return code + [""]

# FIAT_NEW code (compute index function)
# def idx(p,q):
#    return (p+q)*(p+q+1)/2 + q
def _idx(p, q):
    return "((%s) + (%s))*((%s) + (%s) + 1) / 2 + (%s)" % (p, q, p, q, q)

# FIAT_NEW code (helper variables)
# def jrc( a , b , n ):
#    an = float( ( 2*n+1+a+b)*(2*n+2+a+b)) \
#        / float( 2*(n+1)*(n+1+a+b))
#    bn = float( (a*a-b*b) * (2*n+1+a+b) ) \
#        / float( 2*(n+1)*(2*n+a+b)*(n+1+a+b) )
#    cn = float( (n+a)*(n+b)*(2*n+2+a+b)  ) \
#        / float( (n+1)*(n+1+a+b)*(2*n+a+b) )
#    return an,bn,cn
def _jrc(a, b, n):

    an = "( ( 2*(%s)+1+(%s)+(%s))*(2*(%s)+2+(%s)+(%s)) ) / ( 2*((%s)+1)*((%s)+1+(%s)+(%s)))" % (n,a,b,n,a,b)
    bn = "( ((%s)*(%s)-(%s)*(%s)) * (2*(%s)+1+(%s)+(%s))  ) / ( 2*((%s)+1)*(2*(%s)+(%s)+(%s))*((%s)+1+(%s)+(%s)) )" % (a,a,b,b,n,a,b)
    cn = "( ((%s)+(%s))*((%s)+(%s))*(2*(%s)+2+(%s)+(%s))  ) /  ( ((%s)+1)*((%s)+1+(%s)+(%s))*(2*(%s)+(%s)+(%s)) )"  % (n,a,n,b,n,a,b)

    return (an, bn, cn)


def _compute_basisvalues(data, Indent, format):
    """From FIAT_NEW.expansions."""

    code = []

    # Prefetch formats to speed up code generation
    format_add         = format["add"]
    format_multiply         = format["multiply"]
    format_subtract         = format["subtract"]
    format_division         = format["division"]
    format_grouping         = format["grouping"]
    format_sqrt      = format["sqrt"]

#    format_secondary_index  = format["secondary index"]
#    format_psitilde_a       = format["psitilde_a"]
#    format_psitilde_bs      = format["psitilde_bs"]
#    format_psitilde_cs      = format["psitilde_cs"]
#    format_scalings         = format["scalings"]
    format_y                = format["y coordinate"]
#    format_z                = format["z coordinate"]
#    format_const_float      = format["const float declaration"]
    format_float_decl      = format["float declaration"]
    format_basis_table       = format["basisvalues table"]
    format_basisvalue       = format["basisvalues"]
    format_array_access      = format["array access"]
    format_float      = format["floating point"]
    format_uint      = format["uint declaration"]
    format_free_indices      = format["free secondary indices"]
    format_r = format_free_indices[0]
    format_s = format_free_indices[1]

    code += [Indent.indent(format["comment"]("Compute basisvalues"))]

    # Create zero array for basisvalues
    # Get number of members of the expansion set
    # FIXME: KBO: Might be able to get this more elegantly, maybe move to
    # data or simply take the second value of the shape of the coefficients
    embedded_degree = data["embedded_degree"]
    num_mem = data["num_expansion_members"]
    name = format_float_decl + format_basis_table + format_array_access(num_mem)
    value = tabulate_vector([0.0]*num_mem, format)
    code += [(name, value)]

    # Get the element cell domain
    # FIXME: KBO: Change this when supporting R^2 in R^3 elements.
    element_cell_domain = data["cell_domain"]

    # 1D
    if (element_cell_domain == "interval"):
        count = 0
#        for i in range(0, data.degree() + 1):

#            factor = math.sqrt(1.0*i + 0.5)
#            symbol = format_psitilde_a + format_secondary_index(i)

#            # Declare variable
#            name = format_const_float + format_basisvalue(count)

#            # Let inner_product handle format of factor
#            value = inner_product([factor],[symbol], format)

#            code += [(Indent.indent(name), value)]
#            count += 1
#        if (count != poly_dim):
#            error("The number of basis values must be the same as the polynomium dimension of the base")

    # 2D
    elif (element_cell_domain == "triangle"):
        # FIAT_NEW.expansions.TriangleExpansionSet
        # FIXME: KBO: Move common stuff to general functions

        count = 0
        # The initial value basisvalue 0 is always 1.0
        # FIAT_NEW code
        # for ii in range( results.shape[1] ):
        #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
        code += [(format_basisvalue(0), format_float(1.0))]

        # Only continue if the embedded degree is larger than zero
        if embedded_degree > 0:
            # Declare helper variables
            idx0 = "idx0"
            idx1 = "idx1"
            idx2 = "idx2"
            an = "an"
            bn = "bn"
            cn = "cn"
            f1 = "f1"
            f2 = "f2"
            f3 = "f3"
            code += [(format_uint + idx0, 0)]
            code += [(format_uint + idx1, 0)]
            code += [(format_uint + idx2, 0)]
            code += [(format_float_decl + an, format_float(0))]
            code += [(format_float_decl + bn, format_float(0))]
            code += [(format_float_decl + cn, format_float(0))]

            # FIXME: KBO: Use format to compute
            code += [(format_float_decl + f1, "(1.0+2*x+y)/2.0")]
            code += [(format_float_decl + f2, "(1.0 - y) / 2.0")]
            code += [(format_float_decl + f3, "f2*f2")]

            # The initial value of basisfunction 1 is equal to f1
            # FIAT_NEW code
            # results[idx(1,0),:] = f1
            code += [(format_basisvalue(1), f1)]

            # Only active is embedded_degree > 1
            if embedded_degree > 1:
                # FIAT_NEW code
                # for p in range(1,n):
                #    a = (2.0*p+1)/(1.0+p)
                #    b = p / (p+1.0)
                #    results[idx(p+1,0)] = a * f1 * results[idx(p,0),:] \
                #        - p/(1.0+p) * f3 *results[idx(p-1,0),:]
                lines = []
                loop_vars = [(format_r, 1, embedded_degree)]
                lines.append((idx0, _idx("%s + 1" % format_r, 0)))
                lines.append((idx1, _idx(format_r, 0)))
                lines.append((idx2, _idx("%s - 1" % format_r, 0)))
                lines.append((an, "(2.0*%s+1)/(1.0+%s)" % (format_r, format_r)))
                lines.append((bn, "%s/(%s + 1.0)" % (format_r, format_r)))
                fac0 = format_multiply([an, f1, format_basisvalue(idx1)])
                fac1 = format_multiply([format_r + format_division + \
                                        format_grouping(format_add(["1.0", format_r])),\
                                        f3, format_basisvalue(idx2)])
                lines.append((format_basisvalue(idx0), format_subtract([fac0, fac1])))
                code += generate_loop(lines, loop_vars, Indent, format)

                # FIAT_NEW code
                # for p in range(n-1):
                #    for q in range(1,n-p):
                #        (a1,a2,a3) = jrc(2*p+1,0,q)
                #        results[idx(p,q+1),:] \
                #            = ( a1 * y + a2 ) * results[idx(p,q)] \
                #            - a3 * results[idx(p,q-1)]
                lines = []
                loop_vars = [(format_r, 0, embedded_degree - 1),\
                             (format_s, 1, format_subtract([format_float(embedded_degree), format_r]))]
                lines.append((idx0, _idx(format_r, format_add([format_s, format_float(1.0)]))))
                lines.append((idx1, _idx(format_r, format_s)))
                lines.append((idx2, _idx(format_r, format_subtract([format_s, format_float(1.0)]))))
                jrc = _jrc(format_add([format_multiply([format_float(2.0), format_r]), format_float(1.0)]),\
                           format_float(0), format_s)
                lines.append((an, jrc[0]))
                lines.append((bn, jrc[1]))
                lines.append((cn, jrc[2]))
                fac0 = format_multiply([format_grouping(format_add([ \
                        format_multiply([an, format_y]), bn])), format_basisvalue(idx1)])
                fac1 = format_multiply([cn, format_basisvalue(idx2)])
                lines.append((format_basisvalue(idx0), format_subtract([fac0, fac1])))
                code += generate_loop(lines, loop_vars, Indent, format)


            # FIAT_NEW code
            # for p in range(n):
            #    results[idx(p,1),:] = 0.5 * (1+2.0*p+(3.0+2.0*p)*y) \
            #        * results[idx(p,0)]
            lines = []
            loop_vars = [(format_r, 0, embedded_degree)]
            lines.append((idx0, _idx(format_r, 1)))
            lines.append((idx1, _idx(format_r, 0)))
            fac0 = format_multiply([format_float(2.0), format_r])
            fac1 = format_multiply([format_grouping(format_add([format_float(3.0), fac0])), format_y])
            fac2 = format_grouping(format_add([format_float(1.0), fac0, fac1]))
            fac3 = format_multiply([format_float(0.5), fac2, format_basisvalue(idx1)])
            lines.append((format_basisvalue(idx0), fac3))
            code += generate_loop(lines, loop_vars, Indent, format)

            # FIAT_NEW code
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        results[idx(p,q),:] *= math.sqrt((p+0.5)*(p+q+1.0))
            lines = []
            loop_vars = [(format_r, 0, embedded_degree + 1), \
                         (format_s, 0, format_subtract([format_float(embedded_degree + 1), format_r]))]
            lines.append((idx0, _idx(format_r, format_s)))
            fac0 = format_grouping(format_add([format_r, format_float(0.5)]))
            fac1 = format_grouping(format_add([format_r, format_s, format_float(1.0)]))
            fac3 = format_multiply([format_basisvalue(idx0), format_sqrt(format_multiply([fac0, fac1]))])
            lines.append((format_basisvalue(idx0), fac3))
            code += generate_loop(lines, loop_vars, Indent, format)

    # 3D
    elif (element_cell_domain == "tetrahedron"):
        count = 0
#        for k in range(0, data.degree()+1):  # loop over degree
#            for i in range(0, k+1):
#                for j in range(0, k - i + 1):
#                    ii = k-i-j
#                    jj = j
#                    kk = i
#                    factor = math.sqrt( (ii+0.5) * (ii+jj+1.0) * (ii+jj+kk+1.5) )

#                    symbol = format_multiply([format_psitilde_a + format_secondary_index(ii),\
#                             format_scalings(format_y, ii), format_psitilde_bs(ii) + format_secondary_index(jj),\
#                             format_scalings(format_z, (ii+jj)),\
#                             format_psitilde_cs(ii, jj) + format_secondary_index(kk)])

#                    # Declare variable
#                    name = format_const_float + format_basisvalue(count)

#                    # Let inner_product handle format of factor
#                    value = inner_product([factor],[symbol], format)

#                    code += [(Indent.indent(name), value)]
#                    count += 1
#        if (count != poly_dim):
#            error("The number of basis values must be the same as the polynomium dimension of the base")
    else:
        error("Cannot compute basis values for shape: %d" % elemet_cell_domain)

    return code + [""]

#def _eval_jacobi_batch_scalar(a, b, n, variables, format):
#    """Implementation of FIAT function eval_jacobi_batch(a,b,n,xs) from jacobi.py"""

#    # Prefetch formats to speed up code generation
#    format_secondary_index  = format["secondary index"]
#    format_float            = format["floating point"]
#    format_epsilon          = format["epsilon"]

#    # Format variables
#    access = lambda i: variables[1] + format_secondary_index(i)
#    coord = variables[0]

#    if n == 0:
#        return format["floating point"](1.0)
#    if n == 1:
#        # Results for entry 1, of type (a + b * coordinate) (coordinate = x, y or z)
#        res0 = 0.5 * (a - b)
#        res1 = 0.5 * ( a + b + 2.0 )

#        val1 = inner_product([res1], [coord], format)

#        if (abs(res0) > format_epsilon): # Only include if the value is not zero
#            val0 = format_float(res0)
#            if val1:
#                if res0 > 0:
#                    return format["add"]([val1, format_float(res0)])
#                else:
#                    return format["subtract"]([val1, format_float(-res0)])
#            else:
#                return val0
#        else:
#            return val1

#    else:
#        apb = a + b
#        # Compute remaining entries, of type (a + b * coordinate) * psitilde[n-1] - c * psitilde[n-2])
#        a1 = 2.0 * n * ( n + apb ) * ( 2.0 * n + apb - 2.0 )
#        a2 = ( 2.0 * n + apb - 1.0 ) * ( a * a - b * b )
#        a3 = ( 2.0 * n + apb - 2.0 ) * ( 2.0 * n + apb - 1.0 ) * ( 2.0 * n + apb )
#        a4 = 2.0 * ( n + a - 1.0 ) * ( n + b - 1.0 ) * ( 2.0 * n + apb )
#        a2 = a2 / a1
#        a3 = a3 / a1
#        # Note:  we subtract the value of a4!
#        a4 = -a4 / a1

#        float_numbers = [a2, a3, a4]
#        symbols = [access(n-1), format["multiply"]([coord, access(n-1)]), access(n-2)]

#        return inner_product(float_numbers, symbols, format)
