"""Code generation for evaluation of finite element basis values. This module generates
code which is more or less a C++ representation of FIAT code."""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-12-14"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-18

# Python modules
import math
import numpy

# FFC modules
from log import error, debug_code
from removeunused import remove_unused
from codegenerators_utils import IndentControl
from codegenerators_utils import inner_product, tabulate_matrix
from quadrature.quadraturegenerator_utils import generate_loop

# Temporary import
from cpp import format_old as format

def _evaluate_basis(fiat_element):
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

    code = []

    Indent = IndentControl()

    # Get coordinates and generate map
    # FIXME: KBO: If this remains simple, inline function here.
    code += _generate_map(fiat_element, Indent, format)

    # Check if we have just one element

    if (fiat_element.num_sub_elements() == 1):

        # Reset values, change for tensor valued elements
        # FIXME: KBO: If this remains simple, inline function here.
        rank = len(fiat_element.value_shape())
        code += _reset_values(rank, False, Indent, format)

        # Map degree of freedom to local degree of freedom for current element
        code += _map_dof(0, Indent, format)

        # Generate element code
        code += _generate_element_code(fiat_element, 0, False, Indent, format)

    # If the element is of type VectorElement or MixedElement
    else:
        # FIXME: KBO: Only support scalar elements for now
        error("Only scalar elements are supported.")
#        # Reset values, change for tensor valued elements
        # FIXME: KBO: If this remains simple, inline function here.
#        code += _reset_values(fiat_element.value_dimension(0), True, Indent, format)

#        # Generate element code, for all sub-elements
#        code += _mixed_elements(fiat_element, Indent, format)

    lines = format["generate body"](code)

    # Remove unused variables (from transformations and mappings) in code.
    code = remove_unused(lines)
    return code

def _generate_map(fiat_element, Indent, format):
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
#    code += [Indent.indent(format["coordinate map FIAT"](fiat_element.cell_domain()))] + [""]



    # Get coordinates and map to the FIAT reference element from codesnippets.py
    code += [Indent.indent(format["coordinate map FIAT"](cell_domain))] + [""]

    # FIXME: Verify that we don't need to apply additional transformations once we're on the reference element
#    if (fiat_element.cell_domain() == "interval"):

#        # Map coordinates to the reference interval
#        code += [Indent.indent(format_comment("Map coordinates to the reference interval"))]

#        # Code snippet reproduced from FIAT: reference.py: eta_line(xi)
#        code += [Indent.indent(format["snippet eta_interval"])]

#    elif (fiat_element.cell_domain() == "triangle"):

#        # Map coordinates to the reference square
#        code += [Indent.indent(format_comment("Map coordinates to the reference square"))]

#        # Code snippet reproduced from FIAT: reference.py: eta_triangle(xi)
#        code += [Indent.indent(format["snippet eta_triangle"]) %(format_floating_point(format_epsilon))]

#    elif (fiat_element.cell_domain() == "tetrahedron"):

#        # Map coordinates to the reference cube
#        code += [Indent.indent(format_comment("Map coordinates to the reference cube"))]

#        # Code snippet reproduced from FIAT: reference.py: eta_tetrahedron(xi)
#        code += [Indent.indent(format["snippet eta_tetrahedron"]) %(format_floating_point(format_epsilon),\
#                       format_floating_point(format_epsilon))]
#    else:
#        error("Cannot generate map for shape: %d" %(fiat_element.cell_domain()))

    return code

def _reset_values(num_components, vector, Indent, format):
    "Reset all components of the 'values' array as it is a pointer to an array."

    code = []

    if (vector or num_components != 1):
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

def _mixed_elements(fiat_element, Indent, format):
    "Generate code for each sub-element in the event of vector valued elements or mixed elements"

    code = []

    # Prefetch formats to speed up code generation
    format_block_begin = format["block begin"]
    format_block_end = format["block end"]
    format_dof_map_if = format["dof map if"]

    # Extract basis elements, and determine number of elements
    elements = fiat_element.extract_elements()
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

def _generate_element_code(fiat_element, sum_value_dim, vector, Indent, format):
    "Generate code for a single basis element as the dot product of coefficients and basisvalues."

    code = []

#    # Generate basisvalues
#    code += _generate_basisvalues(fiat_element, Indent, format)

    # Tabulate coefficients
    code += _tabulate_coefficients(fiat_element, Indent, format)

    # Compute the value of the basisfunction as the dot product of the coefficients
    # and basisvalues
    code += _compute_values(fiat_element, sum_value_dim, vector, Indent, format)

    return code

#def _generate_basisvalues(fiat_element, Indent, format):
#    "Generate code to compute basis values"

#    code = []

#    # Compute scaling of y and z 1/2(1-y)^n and 1/2(1-z)^n
#    code += _compute_scaling(fiat_element, Indent, format)

#    # Compute auxilliary functions
#    if fiat_element.cell_domain() == "interval":
#        code += _compute_psitilde_a(fiat_element, Indent, format)
#    elif fiat_element.cell_domain() == "triangle":
#        code += _compute_psitilde_a(fiat_element, Indent, format)
#        code += _compute_psitilde_b(fiat_element, Indent, format)
#    elif fiat_element.cell_domain() == "tetrahedron":
#        code += _compute_psitilde_a(fiat_element, Indent, format)
#        code += _compute_psitilde_b(fiat_element, Indent, format)
#        code += _compute_psitilde_c(fiat_element, Indent, format)
#    else:
#        error("Cannot compute auxilliary functions for shape: %d" %(fiat_element.cell_domain()))

#    # Compute the basisvalues
#    code += _compute_basisvalues(fiat_element, Indent, format)

#    return code

def _tabulate_coefficients(fiat_element, Indent, format):
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
    coefficients = fiat_element.get_coeffs()

    # Scalar valued basis element [Lagrange, Discontinuous Lagrange, Crouzeix-Raviart]
    rank = len(fiat_element.value_shape())
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

def _compute_values(fiat_element, sum_value_dim, vector, Indent, format):
    """This function computes the value of the basisfunction as the dot product of the
    coefficients and basisvalues """

    code = []

    # Prefetch formats to speed up code generation
    format_values           = format["argument values"]
    format_array_access     = format["array access"]
    format_matrix_access     = format["matrix access"]
#    format_add              = format["add"]
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
    shape = fiat_element.value_shape()
    if len(shape) > 0:
        num_components = shape[0]
    else:
        num_components = 1

    # Check which transform we should use to map the basis functions
    # FIXME: KBO: Only support affine mapping for now
#    mapping = fiat_element.mapping()
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

#            # Use Piola transform to map basisfunctions back to physical fiat_element if needed
#            if mapping == CONTRAVARIANT_PIOLA:
#                code.insert(i+1,(Indent.indent(format_tmp(0, i)), value))
#                basis_col = [format_tmp_access(0, j) for j in range(fiat_element.cell().topological_dimension())]
#                jacobian_row = [format["transform"]("J", j, i, None) for j in range(fiat_element.cell().topological_dimension())]
#                inner = [format_mult([jacobian_row[j], basis_col[j]]) for j in range(fiat_element.cell().topological_dimension())]
#                sum = format_group(format_add(inner))
#                value = format_mult([format_inv(format_det(None)), sum])
#            elif mapping == COVARIANT_PIOLA:
#                code.insert(i+1,(Indent.indent(format_tmp(0, i)), value))
#                basis_col = [format_tmp_access(0, j) for j in range(fiat_element.cell().topological_dimension())]
#                inverse_jacobian_column = [format["transform"]("JINV", j, i, None) for j in range(fiat_element.cell().topological_dimension())]
#                inner = [format_mult([inverse_jacobian_column[j], basis_col[j]]) for j in range(fiat_element.cell().topological_dimension())]
#                sum = format_group(format_add(inner))
#                value = format_mult([sum])

#            code += [(Indent.indent(name), value)]
    else:
        # Get number of members of the expansion set
        # FIXME: KBO: Might be able to get this more elegantly, maybe move to
        # fiat_element or simply take the second value of the shape of the coefficients
        num_mem = fiat_element.get_nodal_basis().get_expansion_set().get_num_members(fiat_element.get_order())
        loop_vars = [(format_j, 0, num_mem)]

        name = format_pointer + format_values

        value = format_multiply([format_coefficients(0) + format_matrix_access(format_dof, format_j),\
                format_basisvalues(0) + format_array_access(format_j)])

        line = [format["add equal"](name, value)]

        code += generate_loop(line, loop_vars, Indent, format)

    return code

#def _compute_scaling(fiat_element, Indent, format):
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
#    degree = fiat_element.degree()

#    # Get the element cell domain
#    element_cell_domain = fiat_element.cell_domain()

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

#def _compute_psitilde_a(fiat_element, Indent, format):
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
#    degree = fiat_element.degree()

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

#def _compute_psitilde_b(fiat_element, Indent, format):
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
#    degree = fiat_element.degree()

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

#def _compute_psitilde_c(fiat_element, Indent, format):
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
#    degree = fiat_element.degree()

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

#def _compute_basisvalues(fiat_element, Indent, format):
#    """This function is an implementation of the loops inside the FIAT functions
#    tabulate_phis_triangle( n , xs ) and tabulate_phis_tetrahedron( n , xs ) in
#    expansions.py. It computes the basis values from all the previously tabulated variables."""

#    code = []

#    # Prefetch formats to speed up code generation
#    format_multiply         = format["multiply"]
#    format_secondary_index  = format["secondary index"]
#    format_psitilde_a       = format["psitilde_a"]
#    format_psitilde_bs      = format["psitilde_bs"]
#    format_psitilde_cs      = format["psitilde_cs"]
#    format_scalings         = format["scalings"]
#    format_y                = format["y coordinate"]
#    format_z                = format["z coordinate"]
#    format_const_float      = format["const float declaration"]
#    format_basisvalue       = format["basisvalue"]

#    code += [Indent.indent(format["comment"]("Compute basisvalues"))]

#    # Get polynomial dimension of base
#    poly_dim = fiat_element.space_dimension()

#    # Get the element cell domain
#    element_cell_domain = fiat_element.cell_domain()

#    # 1D
#    if (element_cell_domain == "interval"):
#        count = 0
#        for i in range(0, fiat_element.degree() + 1):

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

#    # 2D
#    elif (element_cell_domain == "triangle"):
#        count = 0
#        for k in range(0,fiat_element.degree() + 1):
#            for i in range(0,k + 1):
#                ii = k-i
#                jj = i
#                factor = math.sqrt( (ii+0.5)*(ii+jj+1.0) )

#                symbol = format_multiply([format_psitilde_a + format_secondary_index(ii),\
#                         format_scalings(format_y, ii), format_psitilde_bs(ii) + format_secondary_index(jj)])

#                # Declare variable
#                name = format_const_float + format_basisvalue(count)

#                # Let inner_product handle format of factor
#                value = inner_product([factor],[symbol], format)

#                code += [(Indent.indent(name), value)]
#                count += 1
#        if (count != poly_dim):
#            error("The number of basis values must be the same as the polynomium dimension of the base")

#    # 3D
#    elif (element_cell_domain == "tetrahedron"):
#        count = 0
#        for k in range(0, fiat_element.degree()+1):  # loop over degree
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
#    else:
#        error("Cannot compute basis values for shape: %d" % elemet_cell_domain)

#    return code + [""]

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
