"""Code generation for evaluation of finite element basis values. This module generates
code which is more or less a C++ representation of the code found in FIAT_NEW."""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-04-04"
__copyright__ = "Copyright (C) 2007-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-08

# Python modules
import math
import numpy

# FFC modules
from ffc.log import error, debug_code, ffc_assert
from ffc.cpp import remove_unused
from ffc.cpp import IndentControl
from ffc.cpp import inner_product, tabulate_matrix, tabulate_vector
from ffc.quadrature.quadraturegenerator_utils import generate_loop
from ffc.quadrature.symbolics import create_float
from ffc.quadrature.symbolics import create_symbol
from ffc.quadrature.symbolics import create_sum
from ffc.quadrature.symbolics import create_product
from ffc.quadrature.symbolics import create_fraction
from ffc.quadrature.symbolics import set_format
from ffc.quadrature.symbolics import CONST

# Temporary import
from cpp import format_old as format

def _evaluate_basis(data_list):
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

    set_format(format)
    # FIXME: KBO: Remove when everyting is working
    if data_list == []:
        return ""

    # Init return code and indent object
    code = []
    Indent = IndentControl()

    # Check if we have just one element
    # Get coordinates and generate map
    element_cell_domain = data_list[0]["cell_domain"]
    # FIXME: KBO: This should already be checked elsewhere.
    ffc_assert(all(element_cell_domain == data["cell_domain"] for data in data_list),\
               "The element cell domain must be the same for all sub elements: " + repr(data_list))
    # FIXME: KBO: If this remains simple, inline function here.
    code += _generate_map(element_cell_domain, Indent, format)

    # Get value shape, should also work for TensorElement, scalar are empty tuples, therefore (1,)
    value_shape = sum(sum(data["value_shape"] or (1,)) for data in data_list)
    # Reset values, change for tensor valued elements
    # FIXME: KBO: If this remains simple, inline function here.
    code += _reset_values(value_shape, Indent, format)

    if len(data_list) == 1:
        data = data_list[0]

        # Map degree of freedom to local degree
        code += _map_dof(0, Indent, format)

        # Generate element code
        code += _generate_element_code(data, 0, False, Indent, format)

    # If the element is of type MixedElement (including Vector- and TensorElement)
    else:
        # FIXME: KBO: Only support scalar elements for now
        # Generate element code, for all sub-elements
        code += _mixed_elements(data_list, Indent, format)

    lines = format["generate body"](code)

    # Remove unused variables (from transformations and mappings) in code.
    code = remove_unused(lines)
    return code

def _generate_map(element_cell_domain, Indent, format):
    """Generates map from physical element to the UFC reference element, and from this element
    to reference square/cube. An affine map is assumed for the first mapping and for the second
    map this function implements the UFC version of the FIAT functions, eta_triangle( xi )
    and eta_tetrahedron( xi ) from reference.py"""

    # Prefetch formats to speed up code generation
    format_comment        = format["comment"]
    format_floating_point = format["floating point"]
    format_epsilon        = format["epsilon"]

    # Get the element cell domain
    # FIXME: KBO: Change this when supporting R^2 in R^3 elements.
#    # Get coordinates and map to the UFC reference element from codesnippets.py
#    code += [Indent.indent(format["coordinate map FIAT"](data.cell_domain()))] + [""]

    # Init return code
    code = []

    # Get coordinates and map to the FIAT reference element from codesnippets.py
    code += [Indent.indent(format["coordinate map FIAT"](element_cell_domain))]

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

def _reset_values(value_shape, Indent, format):
    "Reset all components of the 'values' array as it is a pointer to an array."

    # Init return code
    code = []
    if value_shape != 1:
        # Reset values as it is a pointer
        code += [Indent.indent(format["comment"]("Reset values"))]
        code += [(Indent.indent(format["argument values"] + format["array access"](i)),\
                                format["floating point"](0.0)) for i in range(value_shape)]
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

def _mixed_elements(data_list, Indent, format):
    "Generate code for each sub-element in the event of vector valued elements or mixed elements"

    # Prefetch formats to speed up code generation
    format_block_begin = format["block begin"]
    format_block_end = format["block end"]
    format_dof_map_if = format["dof map if"]

    sum_value_dim = 0
    sum_space_dim = 0

    # Init return code
    code = []

    # Generate code for each element
    for data in data_list:

        # Get value and space dimension (should be tensor ready)
        value_dim = sum(data["value_shape"] or (1,))
        space_dim = data["space_dimension"]

        # Determine if the element has a value, for the given dof
        code += [Indent.indent(format_dof_map_if(sum_space_dim, sum_space_dim + space_dim -1))]
        code += [Indent.indent(format_block_begin)]
        # Increase indentation
        Indent.increase()

        # Generate map from global to local dof
        code += _map_dof(sum_space_dim, Indent, format)

        # Generate code for basis element
        code += _generate_element_code(data, sum_value_dim, True, Indent, format)

        # Decrease indentation, finish block - end element code
        Indent.decrease()
        code += [Indent.indent(format_block_end)] + [""]

        # Increase sum of value dimension, and space dimension
        sum_value_dim += value_dim
        sum_space_dim += space_dim

    return code

def _generate_element_code(data, sum_value_dim, vector, Indent, format):
    "Generate code for a single basis element as the dot product of coefficients and basisvalues."

    # Init return code
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

    # Init return code
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

    # Prefetch formats to speed up code generation
    format_comment            = format["comment"]
    format_table_declaration  = format["table declaration"]
    format_coefficients       = format["coefficients table"]
    format_matrix_access      = format["matrix access"]
    format_const_float        = format["const float declaration"]

    # Get coefficients from basis functions, computed by FIAT at compile time
    coefficients = data["coeffs"]

    # Scalar valued basis element [Lagrange, Discontinuous Lagrange, Crouzeix-Raviart]
    family = data["family"]

    # Scalar elements
    if family in ("Lagrange", "Discontinuous Lagrange", "P0"):
        coefficients = [coefficients]
    # Vector valued basis element [Raviart-Thomas, Brezzi-Douglas-Marini (BDM)]
    elif family in ("Brezzi-Douglas-Marini", "Raviart-Thomas", "Nedelec 1st kind H(curl)"):
        coefficients = numpy.transpose(coefficients, [1,0,2])
    # Tensor and other elements
    else:
        error("This finite element family is currently not supported: %s" % family)

    # Init return code
    code = []

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

    # Prefetch formats to speed up code generation
    format_values           = format["argument values"]
    format_array_access     = format["array access"]
    format_matrix_access     = format["matrix access"]
    format_add              = format["add"]

    format_multiply         = format["multiply"]
    format_coefficients     = format["coefficients table"]
    format_basisvalues      = format["basisvalues table"]
    format_j                = format["first free index"]
    format_dof              = format["local dof"]
#    format_secondary_index  = format["secondary index"]
#    format_basisvalue       = format["basisvalue"]
    format_pointer          = format["pointer"]
    format_det              = format["determinant"]
    format_inv              = format["inverse"]
#    format_add              = format["add"]
    format_mult             = format["multiply"]
    format_group            = format["grouping"]
    format_tmp              = format["tmp declaration"]
    format_tmp_access       = format["tmp access"]

    # Init return code
    code = []

    code += [Indent.indent(format["comment"]("Compute value(s)"))]

    # Get number of components, change for tensor valued elements
    shape = data["value_shape"]
    if shape == ():
        num_components = 1
    elif len(shape) == 1:
        num_components = shape[0]
    else:
        error("Tensor valued elements are not supported yet: %s" % data["family"])

    lines = []
    if (vector or num_components != 1):
        # Loop number of components
        for i in range(num_components):
            name = format_values + format_array_access(i + sum_value_dim)

            value = format_multiply([format_coefficients(i) + format_matrix_access(format_dof, format_j),\
                    format_basisvalues + format_array_access(format_j)])
            lines += [format["add equal"](name, value)]
    else:
        # Generate name and value to create matrix vector multiply
        name = format_pointer + format_values
        value = format_multiply([format_coefficients(0) + format_matrix_access(format_dof, format_j),\
                format_basisvalues + format_array_access(format_j)])
        lines = [format["add equal"](name, value)]

    # Get number of members of the expansion set
    num_mem = data["num_expansion_members"]
    loop_vars = [(format_j, 0, num_mem)]
    code += generate_loop(lines, loop_vars, Indent, format)

    # Apply transformatin if applicable
    mapping = data["mapping"]
    if mapping == "affine":
        pass
    elif mapping == "contravariant piola":
        code += ["", Indent.indent(format["comment"]\
                ("Using contravariant Piola transform to map values back to the physical element"))]
        # Get temporary values before mapping
        code += [(Indent.indent(format_tmp(0, i)),\
                  format_values + format_array_access(i + sum_value_dim)) for i in range(num_components)]

        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [format_tmp_access(0, j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create Jacobian.
            jacobian_row = [format["transform"]("J", j, i, None) for j in range(topological_dimension)]

            # Create inner product and multiply by inverse of Jacobian.
            inner = [format_mult([jacobian_row[j], basis_col[j]]) for j in range(topological_dimension)]
            sum_ = format_group(format_add(inner))
            value = format_mult([format_inv(format_det(None)), sum_])
            name = format_values + format_array_access(i + sum_value_dim)
            code += [(name, value)]
    elif mapping == "covariant piola":
        code += ["", Indent.indent(format["comment"]\
                ("Using covariant Piola transform to map values back to the physical element"))]
        # Get temporary values before mapping
        code += [(Indent.indent(format_tmp(0, i)),\
                  format_values + format_array_access(i + sum_value_dim)) for i in range(num_components)]
        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [format_tmp_access(0, j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create inverse of Jacobian.
            inv_jacobian_row = [format["transform"]("JINV", j, i, None) for j in range(topological_dimension)]

            # Create inner product of basis values and inverse of Jacobian.
            inner = [format_mult([inv_jacobian_row[j], basis_col[j]]) for j in range(topological_dimension)]
            value = format_group(format_add(inner))
            name = format_values + format_array_access(i + sum_value_dim)
            code += [(name, value)]
    else:
        error("Unknown mapping: %s" % mapping)

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

# FIAT_NEW code (compute index function) TriangleExpansionSet
# def idx(p,q):
#    return (p+q)*(p+q+1)/2 + q
def _idx2D(p, q):
    f1 = create_float(1)
    f2 = create_float(2)
    idx = create_fraction(create_product([(p+q).expand(), (p+q+f1).expand()]), f2) + q
    return idx

# FIAT_NEW code (compute index function) TetrahedronExpansionSet
# def idx(p,q,r):
#     return (p+q+r)*(p+q+r+1)*(p+q+r+2)/6 + (q+r)*(q+r+1)/2 + r
def _idx3D(p, q, r):
    f1 = create_float(1)
    f2 = create_float(2)
    f6 = create_float(6)
    fac1 = create_fraction( (p + q + r + f2).expand(), f6)
    fac2 = create_fraction( (q + r + f1).expand(), f2)
    fac3 = create_product([(p + q + r).expand(), (p + q + r + f1).expand(), fac1])
    fac4 = create_product([(q + r).expand(), (q + r + f1).expand(), fac2])
    idx = fac3 + fac4 + r
    return idx

# FIAT_NEW code (helper variables) TriangleExpansionSet and TetrahedronExpansionSet
# def jrc( a , b , n ):
#    an = float( ( 2*n+1+a+b)*(2*n+2+a+b)) \
#        / float( 2*(n+1)*(n+1+a+b))

#    bn = float( (a*a-b*b) * (2*n+1+a+b) ) \
#        / float( 2*(n+1)*(2*n+a+b)*(n+1+a+b) )

#    cn = float( (n+a)*(n+b)*(2*n+2+a+b)  ) \
#        / float( (n+1)*(n+1+a+b)*(2*n+a+b) )
#    return an,bn,cn
def _jrc(a, b, n):
    f1 = create_float(1)
    f2 = create_float(2)
    an_num   = create_product([ f2*n + f1 + a + b, f2*n + f2 + a + b ])
    an_denom = create_product([ f2, n + f1, n + f1 + a + b ])
    an = create_fraction(an_num, an_denom)

    bn_num   = create_product([ a*a - b*b, f2*n + f1 + a + b ])
    bn_denom = create_product([ f2, n + f1, f2*n + a + b, n + f1 + a + b ])
    bn = create_fraction(bn_num, bn_denom)

    cn_num   = create_product([ n + a, n + b, f2*n + f2 + a + b ])
    cn_denom = create_product([ n + f1, n + f1 + a + b, f2*n + a + b ])
    cn = create_fraction(cn_num, cn_denom)

    return (an, bn, cn)

def _compute_basisvalues(data, Indent, format):
    """From FIAT_NEW.expansions."""

    # Prefetch formats to speed up code generation
    format_add          = format["add"]
    format_multiply     = format["multiply"]
    format_subtract     = format["subtract"]
    format_division     = format["division"]
    format_grouping     = format["grouping"]
    format_sqrt         = format["sqrt"]
    format_x            = format["x coordinate"]
    format_y            = format["y coordinate"]
    format_z            = format["z coordinate"]
    format_float_decl   = format["float declaration"]
    format_basis_table  = format["basisvalues table"]
    format_basisvalue   = format["basisvalues"]
    format_array_access = format["array access"]
    format_float        = format["floating point"]
    format_uint         = format["uint declaration"]
    format_free_indices = format["free secondary indices"]
    format_r            = format_free_indices[0]
    format_s            = format_free_indices[1]
    format_t            = format_free_indices[2]

    idx0, idx1, idx2    = [format["evaluate_basis aux index"](i) for i in range(1,4)]
    f1, f2, f3, f4, f5  = [create_symbol(format["evaluate_basis aux factor"](i), CONST) for i in range(1,6)]
    an, bn, cn          = [create_symbol(format["evaluate_basis aux value"](i), CONST) for i in range(3)]

    # Get embedded degree
    embedded_degree = data["embedded_degree"]

    # Create helper symbols
    symbol_p = create_symbol(format_r, CONST)
    symbol_q = create_symbol(format_s, CONST)
    symbol_r = create_symbol(format_t, CONST)
    symbol_x = create_symbol(format_x, CONST)
    symbol_y = create_symbol(format_y, CONST)
    symbol_z = create_symbol(format_z, CONST)
    basis_idx0 = create_symbol(format_basisvalue(idx0), CONST)
    basis_idx1 = create_symbol(format_basisvalue(idx1), CONST)
    basis_idx2 = create_symbol(format_basisvalue(idx2), CONST)
    float_0 = create_float(0)
    float_1 = create_float(1)
    float_2 = create_float(2)
    float_3 = create_float(3)
    float_4 = create_float(4)
    float_n = create_float(embedded_degree)
    float_n1 = create_float(embedded_degree + 1)
    float_nm1 = create_float(embedded_degree - 1)
    float_1_5 = create_float(1.5)
    float_0_5 = create_float(0.5)
    float_0_25 = create_float(0.25)

    # Init return code
    code = []

    # Create zero array for basisvalues
    # Get number of members of the expansion set
    num_mem = data["num_expansion_members"]
    name = format_float_decl + format_basis_table + format_array_access(num_mem)
    value = tabulate_vector([0.0]*num_mem, format)
    code += [Indent.indent(format["comment"]("Array of basisvalues"))]
    code += [(name, value)]

    # Declare helper variables, will be removed if not used.
    code += ["", Indent.indent(format["comment"]("Declare helper variables"))]
    code += [(format_uint + idx0, 0)]
    code += [(format_uint + idx1, 0)]
    code += [(format_uint + idx2, 0)]
    code += [(format_float_decl + str(an), format_float(0))]
    code += [(format_float_decl + str(bn), format_float(0))]
    code += [(format_float_decl + str(cn), format_float(0))]

    # Get the element cell domain
    # FIXME: KBO: Change this when supporting R^2 in R^3 elements.
    element_cell_domain = data["cell_domain"]

    # 1D
    if (element_cell_domain == "interval"):
        # FIAT_NEW.expansions.LineExpansionSet
        # FIAT_NEW code
        # psitilde_as = jacobi.eval_jacobi_batch(0,0,n,ref_pts)
        # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs)
        # The initial value basisvalue 0 is always 1.0
        # FIAT_NEW code
        # for ii in range(result.shape[1]):
        #    result[0,ii] = 1.0 + xs[ii,0] - xs[ii,0]
        code += [(format_basisvalue(0), format_float(1.0))]

        # Only continue if the embedded degree is larger than zero
        if embedded_degree > 0:

            # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs)
            # result[1,:] = 0.5 * ( a - b + ( a + b + 2.0 ) * xsnew )
            # The initial value basisvalue 1 is always x
            code += [(format_basisvalue(1), format_x)]
        
            # Only active is embedded_degree > 1
            if embedded_degree > 1:
                # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs)
                # apb = a + b (equal to 0 because of function arguments)
                # for k in range(2,n+1):
                #    a1 = 2.0 * k * ( k + apb ) * ( 2.0 * k + apb - 2.0 )
                #    a2 = ( 2.0 * k + apb - 1.0 ) * ( a * a - b * b )
                #    a3 = ( 2.0 * k + apb - 2.0 )  \
                #        * ( 2.0 * k + apb - 1.0 ) \
                #        * ( 2.0 * k + apb )
                #    a4 = 2.0 * ( k + a - 1.0 ) * ( k + b - 1.0 ) \
                #        * ( 2.0 * k + apb )
                #    a2 = a2 / a1
                #    a3 = a3 / a1
                #    a4 = a4 / a1
                #    result[k,:] = ( a2 + a3 * xsnew ) * result[k-1,:] \
                #        - a4 * result[k-2,:]
                # Declare helper variables (Note the a2 is always zero and therefore left out)
                code += [(format_float_decl + str(f1), float_0)]
                code += [(format_float_decl + str(f2), float_0)]
                code += [(format_float_decl + str(f3), float_0)]
                lines = []
                loop_vars = [(str(symbol_p), 2, float_n1)]
                # Create names
                basis_k = create_symbol(format_basisvalue(str(symbol_p)), CONST)
                basis_km1 = create_symbol(format_basisvalue(str(symbol_p - float_1)), CONST)
                basis_km2 = create_symbol(format_basisvalue(str(symbol_p - float_2)), CONST)
                # Compute helper variables
                a1 = create_product([float_2, symbol_p, symbol_p, float_2*symbol_p - float_2])
                a3 = create_fraction(create_product([float_2*symbol_p,\
                                                     float_2*symbol_p - float_2,
                                                     float_2*symbol_p - float_1]), a1)
                a4 = create_fraction(create_product([float_4*symbol_p, symbol_p - float_1, symbol_p - float_1]), a1)
                lines.append((str(f1), a1 ))
                lines.append((str(f2), a3 ))
                lines.append((str(f3), a4 ))
                # Compute value
                lines.append((str(basis_k), create_product([f2, symbol_x*basis_km1]) - f3*basis_km2))
                # Create loop (block of lines)
                code += generate_loop(lines, loop_vars, Indent, format)

        # Scale values
        # FIAT_NEW.expansions.LineExpansionSet
        # FIAT_NEW code
        # results = numpy.zeros( ( n+1 , len(pts) ) , type( pts[0][0] ) )
        # for k in range( n + 1 ):
        #    results[k,:] = psitilde_as[k,:] * math.sqrt( k + 0.5 )
        lines = []
        loop_vars = [(str(symbol_p), 0, float_n1)]
        # Create names
        basis_k = create_symbol(format_basisvalue(str(symbol_p)), CONST)
        # Compute value
        fac1 = create_symbol( format_sqrt(str(symbol_p + float_0_5)), CONST )
        lines += [format["times equal"](str(basis_k), str(fac1))]
        # Create loop (block of lines)
        code += generate_loop(lines, loop_vars, Indent, format)
    # 2D
    elif (element_cell_domain == "triangle"):
        # FIAT_NEW.expansions.TriangleExpansionSet

        # Compute helper factors
        # FIAT_NEW code
        # f1 = (1.0+2*x+y)/2.0
        # f2 = (1.0 - y) / 2.0
        # f3 = f2**2
        fac1 = create_fraction(float_1 + float_2*symbol_x + symbol_y, float_2)
        fac2 = create_fraction(float_1 - symbol_y, float_2)
        code += [(format_float_decl + str(f1), fac1)]
        code += [(format_float_decl + str(f2), fac2)]
        code += [(format_float_decl + str(f3), f2*f2)]

        code += ["", Indent.indent(format["comment"]("Compute basisvalues"))]
        # The initial value basisvalue 0 is always 1.0
        # FIAT_NEW code
        # for ii in range( results.shape[1] ):
        #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
        code += [(format_basisvalue(0), format_float(1.0))]

        # Only continue if the embedded degree is larger than zero
        if embedded_degree > 0:

            # The initial value of basisfunction 1 is equal to f1
            # FIAT_NEW code
            # results[idx(1,0),:] = f1
            code += [(format_basisvalue(1), str(f1))]

            # Only active is embedded_degree > 1
            if embedded_degree > 1:
                # FIAT_NEW code
                # for p in range(1,n):
                #    a = (2.0*p+1)/(1.0+p)
                #    b = p / (p+1.0)
                #    results[idx(p+1,0)] = a * f1 * results[idx(p,0),:] \
                #        - p/(1.0+p) * f3 *results[idx(p-1,0),:]
                # FIXME: KBO: Is there an error in FIAT? why is b not used?
                lines = []
                loop_vars = [(str(symbol_p), 1, embedded_degree)]
                # Compute indices
                lines.append((idx0, _idx2D(symbol_p + float_1, float_0)))
                lines.append((idx1, _idx2D(symbol_p, float_0)))
                lines.append((idx2, _idx2D(symbol_p - float_1, float_0)))
                # Compute single helper variable an
                lines.append((str(an), create_fraction(float_2*symbol_p + float_1, symbol_p + float_1)))
                # Compute value
                fac0 = create_product([an, f1, basis_idx1])
                fac1 = create_product([create_fraction(symbol_p, float_1 + symbol_p), f3, basis_idx2])
                lines.append((str(basis_idx0), fac0 - fac1))
                # Create loop (block of lines)
                code += generate_loop(lines, loop_vars, Indent, format)

                # FIAT_NEW code
                # for p in range(n-1):
                #    for q in range(1,n-p):
                #        (a1,a2,a3) = jrc(2*p+1,0,q)
                #        results[idx(p,q+1),:] \
                #            = ( a1 * y + a2 ) * results[idx(p,q)] \
                #            - a3 * results[idx(p,q-1)]
                lines = []
                loop_vars = [(str(symbol_p), 0, embedded_degree - 1),\
                             (str(symbol_q), 1, float_n - symbol_p)]
                # Compute indices
                lines.append((idx0, _idx2D(symbol_p, symbol_q + float_1)))
                lines.append((idx1, _idx2D(symbol_p, symbol_q)))
                lines.append((idx2, _idx2D(symbol_p, symbol_q - float_1)))
                # Comute all helper variables
                jrc = _jrc(float_2*symbol_p + float_1, float_0, symbol_q)
                lines.append((str(an), jrc[0]))
                lines.append((str(bn), jrc[1]))
                lines.append((str(cn), jrc[2]))
                # Compute value
                fac0 = create_product([an * symbol_y + bn, basis_idx1])
                fac1 = cn * basis_idx2
                lines.append((str(basis_idx0), fac0 - fac1))
                # Create loop (block of lines)
                code += generate_loop(lines, loop_vars, Indent, format)

            # FIAT_NEW code
            # for p in range(n):
            #    results[idx(p,1),:] = 0.5 * (1+2.0*p+(3.0+2.0*p)*y) \
            #        * results[idx(p,0)]
            lines = []
            loop_vars = [(str(symbol_p), 0, embedded_degree)]
            # Compute indices
            lines.append((idx0, _idx2D(symbol_p, float_1)))
            lines.append((idx1, _idx2D(symbol_p, float_0)))
            # Compute value
            fac0 = create_product([float_3 + float_2*symbol_p, symbol_y])
            fac1 = create_product([float_0_5, float_1 + float_2*symbol_p + fac0, basis_idx1])
            lines.append((str(basis_idx0), fac1.expand().reduce_ops()))
            # Create loop (block of lines)
            code += generate_loop(lines, loop_vars, Indent, format)

            # FIAT_NEW code
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        results[idx(p,q),:] *= math.sqrt((p+0.5)*(p+q+1.0))
            lines = []
            loop_vars = [(str(symbol_p), 0, embedded_degree + 1), \
                         (str(symbol_q), 0, float_n1 - symbol_p)]
            # Compute indices
            lines.append((idx0, _idx2D(symbol_p, symbol_q)))
            # Compute value
            fac0 = create_product([symbol_p + float_0_5, symbol_p + symbol_q + float_1])
            fac2 = create_symbol( format_sqrt(str(fac0)), CONST )
            lines += [format["times equal"](str(basis_idx0), fac2)]
            # Create loop (block of lines)
            code += generate_loop(lines, loop_vars, Indent, format)

    # 3D
    elif (element_cell_domain == "tetrahedron"):
        # FIAT_NEW.expansions.TetrahedronExpansionSet

        # Compute helper factors
        # FIAT_NEW code
        # factor1 = 0.5 * ( 2.0 + 2.0*x + y + z )
        # factor2 = (0.5*(y+z))**2
        # factor3 = 0.5 * ( 1 + 2.0 * y + z )
        # factor4 = 0.5 * ( 1 - z )
        # factor5 = factor4 ** 2
        fac1 = create_product([float_0_5, float_2 + float_2*symbol_x + symbol_y + symbol_z])
        fac2 = create_product([float_0_25, symbol_y + symbol_z, symbol_y + symbol_z])
        fac3 = create_product([float_0_5, float_1 + float_2*symbol_y + symbol_z])
        fac4 = create_product([float_0_5, float_1 - symbol_z])
        code += [(format_float_decl + str(f1), fac1)]
        code += [(format_float_decl + str(f2), fac2)]
        code += [(format_float_decl + str(f3), fac3)]
        code += [(format_float_decl + str(f4), fac4)]
        code += [(format_float_decl + str(f5), f4*f4)]

        code += ["", Indent.indent(format["comment"]("Compute basisvalues"))]
        # The initial value basisvalue 0 is always 1.0
        # FIAT_NEW code
        # for ii in range( results.shape[1] ):
        #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
        code += [(format_basisvalue(0), format_float(1.0))]

        # Only continue if the embedded degree is larger than zero
        if embedded_degree > 0:

            # The initial value of basisfunction 1 is equal to f1
            # FIAT_NEW code
            # results[idx(1,0),:] = f1
            code += [(format_basisvalue(1), str(f1))]

            # Only active is embedded_degree > 1
            if embedded_degree > 1:

                # FIAT_NEW code
                # for p in range(1,n):
                #    a1 = ( 2.0 * p + 1.0 ) / ( p + 1.0 )
                #    a2 = p / (p + 1.0)
                #    results[idx(p+1,0,0)] = a1 * factor1 * results[idx(p,0,0)] \
                #        -a2 * factor2 * results[ idx(p-1,0,0) ]
                lines = []
                loop_vars = [(str(symbol_p), 1, embedded_degree)]
                # Compute indices
                lines.append((idx0, _idx3D(symbol_p + float_1, float_0, float_0)))
                lines.append((idx1, _idx3D(symbol_p          , float_0, float_0)))
                lines.append((idx2, _idx3D(symbol_p - float_1, float_0, float_0)))
                # Compute value
                fac1 = create_fraction(float_2*symbol_p + float_1, symbol_p + float_1)
                fac2 = create_fraction(symbol_p, symbol_p + float_1)
                fac3 = create_product([fac1, f1, basis_idx1]) - create_product([fac2, f2, basis_idx2]) 
                lines.append((str(basis_idx0), fac3))
                # Create loop (block of lines)
                code += generate_loop(lines, loop_vars, Indent, format)

                # FIAT_NEW code
                # for p in range(0,n-1):
                #    for q in range(1,n-p):
                #        (aq,bq,cq) = jrc(2*p+1,0,q)
                #        qmcoeff = aq * factor3 + bq * factor4
                #        qm1coeff = cq * factor5
                #        results[idx(p,q+1,0)] = qmcoeff * results[idx(p,q,0)] \
                #            - qm1coeff * results[idx(p,q-1,0)]
                lines = []
                loop_vars = [(str(symbol_p), 0, embedded_degree - 1),\
                             (str(symbol_q), 1, float_n - symbol_p)]
                # Compute indices
                lines.append((idx0, _idx3D(symbol_p, symbol_q + float_1, float_0)))
                lines.append((idx1, _idx3D(symbol_p, symbol_q          , float_0)))
                lines.append((idx2, _idx3D(symbol_p, symbol_q - float_1, float_0)))
                # Comute all helper variables
                jrc = _jrc(float_2*symbol_p + float_1, float_0, symbol_q)
                lines.append((str(an), jrc[0]))
                lines.append((str(bn), jrc[1]))
                lines.append((str(cn), jrc[2]))
                # Compute value
                fac1 = create_product([an*f3 + bn*f4, basis_idx1]) - cn*f5*basis_idx2
                lines.append((str(basis_idx0), fac1))
                # Create loop (block of lines)
                code += generate_loop(lines, loop_vars, Indent, format)

                # FIAT_NEW code
                # general r by recurrence
                # for p in range(n-1):
                #     for q in range(0,n-p-1):
                #         for r in range(1,n-p-q):
                #             ar,br,cr = jrc(2*p+2*q+2,0,r)
                #             results[idx(p,q,r+1)] = \
                #                         (ar * z + br) * results[idx(p,q,r) ] \
                #                         - cr * results[idx(p,q,r-1) ]
                lines = []
                loop_vars = [(str(symbol_p), 0, embedded_degree - 1),\
                             (str(symbol_q), 0, float_nm1 - symbol_p),\
                             (str(symbol_r), 1, float_n - symbol_p - symbol_q)]
                # Compute indices
                lines.append((idx0, _idx3D(symbol_p, symbol_q, symbol_r + float_1)))
                lines.append((idx1, _idx3D(symbol_p, symbol_q, symbol_r)))
                lines.append((idx2, _idx3D(symbol_p, symbol_q, symbol_r - float_1)))
                # Comute all helper variables
                jrc = _jrc(float_2*symbol_p + float_2*symbol_q, float_0, symbol_r)
                lines.append((str(an), jrc[0]))
                lines.append((str(bn), jrc[1]))
                lines.append((str(cn), jrc[2]))
                # Compute value
                fac1 = create_product([an*symbol_z + bn, basis_idx1]) - cn*basis_idx2
                lines.append((str(basis_idx0), fac1))
                # Create loop (block of lines)
                code += generate_loop(lines, loop_vars, Indent, format)

            # FIAT_NEW code
            # q = 1
            # for p in range(0,n):
            #    results[idx(p,1,0)] = results[idx(p,0,0)] \
            #        * ( p * (1.0 + y) + ( 2.0 + 3.0 * y + z ) / 2 )
            lines = []
            loop_vars = [(str(symbol_p), 0, embedded_degree)]
            # Compute indices
            lines.append((idx0, _idx3D(symbol_p, float_1, float_0)))
            lines.append((idx1, _idx3D(symbol_p, float_0, float_0)))
            # Compute value
            fac1 = create_fraction(float_2 + float_3*symbol_y + symbol_z, float_2)
            fac2 = create_product([symbol_p, float_1 + symbol_y])
            fac3 = create_product([basis_idx1, fac2 + fac1])
            lines.append((str(basis_idx0), fac3))
            # Create loop (block of lines)
            code += generate_loop(lines, loop_vars, Indent, format)
                
            # FIAT_NEW code
            # now handle r=1
            # for p in range(n):
            #    for q in range(n-p):
            #        results[idx(p,q,1)] = results[idx(p,q,0)] \
            #            * ( 1.0 + p + q + ( 2.0 + q + p ) * z )
            lines = []
            loop_vars = [(str(symbol_p), 0, embedded_degree),\
                         (str(symbol_q), 0, float_n - symbol_p)]
            # Compute indices
            lines.append((idx0, _idx3D(symbol_p, symbol_q, float_1)))
            lines.append((idx1, _idx3D(symbol_p, symbol_q, float_0)))
            # Compute value
            fac1 = create_product([float_2 + symbol_p + symbol_q, symbol_z])
            fac2 = create_product([basis_idx1, float_1 + symbol_p + symbol_q + fac1])
            lines.append((str(basis_idx0), fac2))
            # Create loop (block of lines)
            code += generate_loop(lines, loop_vars, Indent, format)

            # FIAT_NEW code
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        for r in range(n-p-q+1):
            #            results[idx(p,q,r)] *= math.sqrt((p+0.5)*(p+q+1.0)*(p+q+r+1.5))
            lines = []
            loop_vars = [(str(symbol_p), 0, float_n1),\
                         (str(symbol_q), 0, float_n1 - symbol_p),\
                         (str(symbol_r), 0, float_n1 - symbol_p - symbol_q)]
            # Compute indices
            lines.append((idx0, _idx3D(symbol_p, symbol_q, symbol_r)))
            # Compute value
            fac0 = create_product([symbol_p + float_0_5,\
                                   symbol_p + symbol_q + float_1,\
                                   symbol_p + symbol_q + symbol_r + float_1_5])
            fac2 = create_symbol( format_sqrt(str(fac0)), CONST )
            lines += [format["times equal"](str(basis_idx0), fac2)]
            # Create loop (block of lines)
            code += generate_loop(lines, loop_vars, Indent, format)
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
