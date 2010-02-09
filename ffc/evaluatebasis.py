"""Code generation for evaluation of finite element basis values. This module generates
code which is more or less a C++ representation of the code found in FIAT."""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2007-04-04"
__copyright__ = "Copyright (C) 2007-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-02-04

# Python modules
import math
import numpy

# FFC modules
from ffc.log import error
from ffc.cpp import remove_unused, indent, format
from ffc.quadrature.symbolics import create_float, create_float, create_symbol,\
                                     create_product, create_sum, create_fraction, CONST

def _evaluate_basis_all(data_list):
    """Like evaluate_basis, but return the values of all basis functions (dofs)."""

    if isinstance(data_list, str):
        return format["exception"]("evaluate_basis_all: %s" % data_list)

    # Prefetch formats.
    f_assign    = format["assign"]
    f_component = format["component"]
    f_comment   = format["comment"]
    f_loop      = format["generate loop"]
    f_r, f_s    = format["free indices"][:2]
    f_tensor    = format["tabulate tensor"]
    f_values    = format["argument values"]
    f_basis     = format["call basis"]
    f_dof_vals  = format["dof values"]
    f_double    = format["float declaration"]
    f_float     = format["floating point"]
    f_decl      = format["declaration"]
    f_ref_var   = format["reference variable"]

    # Initialise return code.
    code = []

    # FIXME: KBO: Figure out how the return format should be, either:
    # [N0[0], N0[1], N1[0], N1[1], ...]
    # or
    # [N0[0], N1[0], ..., N0[1], N1[1], ...]
    # for vector (tensor elements), currently returning option 1.

    # FIXME: KBO: For now, just call evaluate_basis and map values accordingly,
    # this will keep the amount of code at a minimum. If it turns out that speed
    # is an issue (overhead from calling evaluate_basis), we can easily generate
    # all the code

    # Get total value shape and space dimension for entire element (possibly mixed).
    value_shape = sum(sum(data["value_shape"] or (1,)) for data in data_list)
    space_dimension = sum(data["space_dimension"] for data in data_list)

    # Special case where space dimension is one (constant elements).
    if space_dimension == 1:
        code += [f_comment("Element is constant, calling evaluate_basis.")]
        code += [f_basis(format["int"](0), f_values)]
        return "\n".join(code)

    # Declare helper value to hold single dof values.
    code += [f_comment("Helper variable to hold values of a single dof.")]
    if value_shape == 1:
        code += [f_decl(f_double, f_dof_vals, f_float(0.0))]
    else:
        code += [f_decl(f_double, f_component(f_dof_vals, value_shape), f_tensor([0.0]*value_shape))]

    # Create loop over dofs that calls evaluate_basis for a single dof and
    # inserts the values into the global array.
    code += ["", f_comment("Loop dofs and call evaluate_basis.")]
    lines_r = []
    loop_vars_r = [(f_r, 0, space_dimension)]

    if value_shape == 1:
        lines_r += [f_basis(f_r, f_ref_var(f_dof_vals))]
    else:
        lines_r += [f_basis(f_r, f_dof_vals)]

    if value_shape ==  1:
        lines_r += [f_assign(f_component(f_values, f_r), f_dof_vals)]
    else:
        index = format["matrix index"](f_r, f_s, value_shape)
        lines_s = [f_assign(f_component(f_values, index), f_component(f_dof_vals, f_s))]
        lines_r += f_loop(lines_s, [(f_s, 0, value_shape)])

    code += f_loop(lines_r, loop_vars_r)

    # Generate bode (no need to remove unused).
    return "\n".join(code)


# From FIAT_NEW.polynomial_set.tabulate()
def _evaluate_basis(data_list):
    """Generate run time code to evaluate an element basisfunction at an
    arbitrary point. The value(s) of the basisfunction is/are
    computed as in FIAT as the dot product of the coefficients (computed at compile time)
    and basisvalues which are dependent on the coordinate and thus have to be computed at
    run time.

    The function should work for all elements supported by FIAT, but it remains
    untested for tensor valued elements."""

    if isinstance(data_list, str):
        return format["exception"]("evaluate_basis: %s" % data_list)

    # Prefetch formats.
    f_assign    = format["assign"]
    f_comment   = format["comment"]
    f_values    = format["argument values"]
    f_float     = format["floating point"]
    f_component = format["component"]

    # Initialise return code.
    code = []

    # Get the element cell domain and geometric dimension.
    element_cell_domain = data_list[0]["cell_domain"]
    geometric_dimension = data_list[0]["geometric_dimension"]

    # Get code snippets for Jacobian, Inverse of Jacobian and mapping of
    # coordinates from physical element to the FIAT reference element.
    # FIXME: KBO: Change this when supporting R^2 in R^3 elements.
    code += [format["jacobian and inverse"](geometric_dimension)]
    code += ["", format["fiat coordinate map"](element_cell_domain)]

    # Get value shape and reset values. This should also work for TensorElement,
    # scalar are empty tuples, therefore (1,) in which case value_shape = 1.
    value_shape = sum(sum(data["value_shape"] or (1,)) for data in data_list)
    code += ["", f_comment("Reset values.")]
    if value_shape == 1:
        # Reset values as a pointer.
        code += [f_assign(format["dereference pointer"](f_values), f_float(0.0))]
    else:
        # Reset all values.
        code += [f_assign(f_component(f_values, i), f_float(0.0)) for i in range(value_shape)]

    if len(data_list) == 1:
        data = data_list[0]

        # Map degree of freedom to local degree.
        code += ["", _map_dof(0)]

        # Generate element code.
        code += _generate_element_code(data, 0, False)

    # If the element is of type MixedElement (including Vector- and TensorElement).
    else:
        # Generate element code, for all sub-elements.
        code += _mixed_elements(data_list)

    # Remove unused variables (from transformations and mappings) in code.
    code = remove_unused("\n".join(code))
    return code

def _map_dof(sum_space_dim):
    """This function creates code to map a basis function to a local basis function.
    Example, the following mixed element:

    element = VectorElement("Lagrange", "triangle", 2)

    has the element list, elements = [Lagrange order 2, Lagrange order 2] and 12 dofs (6 each).
    The evaluation of basis function 8 is then mapped to 2 (8-6) for local element no. 2."""

    # In case of only one element or the first element in a series then we don't subtract anything.
    if sum_space_dim == 0:
        code = [format["comment"]("Map degree of freedom to element degree of freedom")]
        code += [format["const uint declaration"](format["local dof"], format["argument basis num"])]
    else:
        code = [format["comment"]("Map degree of freedom to element degree of freedom")]
        code += [format["const uint declaration"](format["local dof"],\
                format["sub"]([format["argument basis num"], "%d" % sum_space_dim]))]

    return "\n".join(code)

def _mixed_elements(data_list):
    "Generate code for each sub-element in the event of mixed elements"

    # Prefetch formats to speed up code generation.
    f_dof_map_if = format["dof map if"]
    f_if         = format["if"]

    sum_value_dim = 0
    sum_space_dim = 0

    # Initialise return code.
    code = []

    # Loop list of data and generate code for each element.
    for data in data_list:

        # Get value and space dimension (should be tensor ready).
        value_dim = sum(data["value_shape"] or (1,))
        space_dim = data["space_dimension"]

        # Generate map from global to local dof.
        element_code = [_map_dof(sum_space_dim)]

        # Generate code for basis element.
        element_code += _generate_element_code(data, sum_value_dim, True)

        # Remove unused code for each sub element and indent code.
        if_code = indent(remove_unused("\n".join(element_code)), 2)

        # Create if statement and add to code.
        code += [f_if(f_dof_map_if(sum_space_dim, sum_space_dim + space_dim - 1), if_code)]

        # Increase sum of value dimension, and space dimension.
        sum_value_dim += value_dim
        sum_space_dim += space_dim

    return code

def _generate_element_code(data, sum_value_dim, vector):
    """Generate code for a single basis element as the dot product of
    coefficients and basisvalues. Then apply transformation if applicable."""

    # Generate basisvalues.
    code = _compute_basisvalues(data)

    # Tabulate coefficients.
    code += _tabulate_coefficients(data)

    # Compute the value of the basisfunction as the dot product of the coefficients
    # and basisvalues and apply transformation.
    code += _compute_values(data, sum_value_dim, vector)

    return code

def _tabulate_coefficients(data):
    """This function tabulates the element coefficients that are generated by FIAT at
    compile time."""

    # Prefetch formats to speed up code generation.
    f_comment       = format["comment"]
    f_table         = format["static const float declaration"]
    f_coefficients  = format["coefficients"]
    f_component     = format["component"]
    f_decl          = format["declaration"]
    f_tensor        = format["tabulate tensor"]
    f_new_line      = format["new line"]

    # Get coefficients from basis functions, computed by FIAT at compile time.
    coefficients = data["coeffs"]

    # Use rank to handle coefficients.
    rank = len(data["value_shape"])
    if rank == 0:
        coefficients = [coefficients]
    # Vector valued basis element [Raviart-Thomas, Brezzi-Douglas-Marini (BDM)].
    elif rank == 1:
        coefficients = numpy.transpose(coefficients, [1,0,2])
    # Tensor and other elements.
    else:
        error("Rank %d elements are currently not supported" % rank)

    # Initialise return code.
    code = [f_comment("Table(s) of coefficients.")]

    # Generate tables for each component.
    for i, coeffs in enumerate(coefficients):
        # Get number of dofs and number of members of the expansion set.
        num_dofs, num_mem = numpy.shape(coeffs)

        # Varable name for coefficients.
        name = f_component(f_coefficients(i), [num_dofs, num_mem])

        # Generate array of values.
        code += [f_decl(f_table, name, f_new_line + f_tensor(coeffs))] + [""]
    return code

def _compute_values(data, sum_value_dim, vector):
    """This function computes the value of the basisfunction as the dot product
    of the coefficients and basisvalues."""

    # Prefetch formats to speed up code generation.
    f_values        = format["argument values"]
    f_component     = format["component"]
    f_comment       = format["comment"]
    f_add           = format["add"]
    f_coefficients  = format["coefficients"]
    f_basisvalues   = format["basisvalues"]
    f_r             = format["free indices"][0]
    f_dof           = format["local dof"]
    f_deref_pointer = format["dereference pointer"]
    f_detJ          = format["det(J)"]
    f_inv           = format["inverse"]
    f_mul           = format["mul"]
    f_iadd          = format["iadd"]
    f_group         = format["grouping"]
    f_tmp_ref       = format["tmp ref value"]
    f_assign        = format["assign"]
    f_loop          = format["generate loop"]
    f_const_float   = format["const float declaration"]
    f_trans         = format["transform"]
    f_inner         = format["inner product"]

    # Initialise return code.
    code = [f_comment("Compute value(s).")]

    # Get number of components, change for tensor valued elements.
    shape = data["value_shape"]
    if shape == ():
        num_components = 1
    elif len(shape) == 1:
        num_components = shape[0]
    else:
        error("Tensor valued elements are not supported yet: %d " % shape)

    lines = []
    if (vector or num_components != 1):
        # Loop number of components.
        for i in range(num_components):
            # Generate name and value to create matrix vector multiply.
            name = f_component(f_values, i + sum_value_dim)
            value = f_mul([f_component(f_coefficients(i), [f_dof, f_r]),\
                    f_component(f_basisvalues, f_r)])
            lines += [f_iadd(name, value)]
    else:
        # Generate name and value to create matrix vector multiply.
        name = f_deref_pointer(f_values)
        value = f_mul([f_component(f_coefficients(0), [f_dof, f_r]),\
                f_component(f_basisvalues, f_r)])
        lines = [f_iadd(name, value)]

    # Get number of members of the expansion set.
    num_mem = data["num_expansion_members"]
    loop_vars = [(f_r, 0, num_mem)]
    code += f_loop(lines, loop_vars)

    # Apply transformation if applicable.
    mapping = data["mapping"]
    if mapping == "affine":
        pass
    elif mapping == "contravariant piola":
        code += ["", f_comment("Using contravariant Piola transform to map values back to the physical element.")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + sum_value_dim))\
                  for i in range(num_components)]
        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create Jacobian.
            jacobian_row = [f_trans("J", i, j, None) for j in range(topological_dimension)]

            # Create inner product and multiply by inverse of Jacobian.
            inner = f_group(f_inner(jacobian_row, basis_col))
            value = f_mul([f_inv(f_detJ(None)), inner])
            name = f_component(f_values, i + sum_value_dim)
            code += [f_assign(name, value)]
    elif mapping == "covariant piola":
        code += ["", f_comment("Using covariant Piola transform to map values back to the physical element.")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + sum_value_dim))\
                  for i in range(num_components)]
        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create inverse of Jacobian.
            inv_jacobian_column = [f_trans("JINV", j, i, None) for j in range(topological_dimension)]

            # Create inner product of basis values and inverse of Jacobian.
            value = f_group(f_inner(inv_jacobian_column, basis_col))
            name = f_component(f_values, i + sum_value_dim)
            code += [f_assign(name, value)]
    else:
        error("Unknown mapping: %s" % mapping)

    return code

# FIAT_NEW code (compute index function) TriangleExpansionSet.
# def idx(p,q):
#    return (p+q)*(p+q+1)/2 + q
def _idx2D(p, q):
    f_add   = format["add"]
    f_group = format["grouping"]
    f_int   = format["int"]
    pq = format["addition"]([str(p), str(q)])
    pq1 = f_group(f_add([pq, f_int(1)]))
    if q == "0":
        return format["div"](format["mul"]([pq, pq1]), f_int(2))
    return f_add([format["div"](format["mul"]([f_group(pq), pq1]), f_int(2)), str(q)])

# FIAT_NEW code (compute index function) TetrahedronExpansionSet.
# def idx(p,q,r):
#     return (p+q+r)*(p+q+r+1)*(p+q+r+2)/6 + (q+r)*(q+r+1)/2 + r
def _idx3D(p, q, r):
    f_add   = format["add"]
    f_group = format["grouping"]
    f_int   = format["int"]
    pqr  = format["addition"]([str(p), str(q), str(r)])
    pqr1 = f_group(f_add([pqr, f_int(1)]))
    pqr2 = f_group(f_add([pqr, f_int(2)]))
    qr   = format["addition"]([str(q), str(r)])
    qr1  = f_group(f_add([qr, f_int(1)]))
    if q == r == "0":
        return format["div"](format["mul"]([pqr, pqr1, pqr2]), f_int(6))
    pqrg = f_group(pqr)
    fac0 = format["div"](format["mul"]([pqrg, pqr1, pqr2]), f_int(6))
    if r == "0":
        return f_add([fac0, format["div"](format["mul"]([qr, qr1]), f_int(2))])
    return f_add([fac0, format["div"](format["mul"]([f_group(qr), qr1]), f_int(2)), str(r)])

# FIAT_NEW code (helper variables) TriangleExpansionSet and TetrahedronExpansionSet.
# def jrc( a , b , n ):
#    an = float( ( 2*n+1+a+b)*(2*n+2+a+b)) \
#        / float( 2*(n+1)*(n+1+a+b))

#    bn = float( (a*a-b*b) * (2*n+1+a+b) ) \
#        / float( 2*(n+1)*(2*n+a+b)*(n+1+a+b) )

#    cn = float( (n+a)*(n+b)*(2*n+2+a+b)  ) \
#        / float( (n+1)*(n+1+a+b)*(2*n+a+b) )
#    return an,bn,cn
def _jrc(a, b, n):
    f1        = create_float(1)
    f2        = create_float(2)
    an_num    = create_product([ f2*n + f1 + a + b, f2*n + f2 + a + b ])
    an_denom  = create_product([ f2, n + f1, n + f1 + a + b ])
    an        = create_fraction(an_num, an_denom)

    bn_num    = create_product([ a*a - b*b, f2*n + f1 + a + b ])
    bn_denom  = create_product([ f2, n + f1, f2*n + a + b, n + f1 + a + b ])
    bn        = create_fraction(bn_num, bn_denom)

    cn_num    = create_product([ n + a, n + b, f2*n + f2 + a + b ])
    cn_denom  = create_product([ n + f1, n + f1 + a + b, f2*n + a + b ])
    cn        = create_fraction(cn_num, cn_denom)
    return (an, bn, cn)

def _compute_basisvalues(data):
    """From FIAT_NEW.expansions."""

    # Prefetch formats to speed up code generation.
    f_comment     = format["comment"]
    f_add         = format["add"]
    f_imul        = format["imul"]
    f_sub         = format["sub"]
    f_group       = format["grouping"]
    f_assign      = format["assign"]
    f_sqrt        = format["sqrt"]
    f_x           = format["x coordinate"]
    f_y           = format["y coordinate"]
    f_z           = format["z coordinate"]
    f_double      = format["float declaration"]
    f_basisvalue  = format["basisvalues"]
    f_component   = format["component"]
    f_float       = format["floating point"]
    f_uint        = format["uint declaration"]
    f_tensor      = format["tabulate tensor"]
    f_loop        = format["generate loop"]
    f_decl        = format["declaration"]
    f_tmp         = format["tmp value"]
    f_int         = format["int"]

    f_r, f_s, f_t = format["free indices"][:3]
    idx0          = f_r + f_r
    idx1          = f_s + f_s
    idx2          = f_t + f_t

    # Create temporary values.
    f1, f2, f3, f4, f5  = [create_symbol(f_tmp(i), CONST) for i in range(0,5)]
    an, bn, cn          = [create_symbol(f_tmp(i), CONST) for i in range(5,8)]

    # Get embedded degree.
    embedded_degree = data["embedded_degree"]

    # Create helper symbols.
    symbol_p    = create_symbol(f_r, CONST)
    symbol_q    = create_symbol(f_s, CONST)
    symbol_r    = create_symbol(f_t, CONST)
    symbol_x    = create_symbol(f_x, CONST)
    symbol_y    = create_symbol(f_y, CONST)
    symbol_z    = create_symbol(f_z, CONST)
    basis_idx0  = create_symbol(f_component(f_basisvalue, idx0), CONST)
    basis_idx1  = create_symbol(f_component(f_basisvalue, idx1), CONST)
    basis_idx2  = create_symbol(f_component(f_basisvalue, idx2), CONST)
    int_0     = f_int(0)
    int_1     = f_int(1)
    int_2    = f_int(2)
    int_n    = f_int(embedded_degree)
    int_n1   = f_int(embedded_degree + 1)
    int_nm1  = f_int(embedded_degree - 1)
    float_0 = create_float(0)
    float_1 = create_float(1)
    float_2 = create_float(2)
    float_3 = create_float(3)
    float_4 = create_float(4)
    float_1_5   = create_float(1.5)
    float_0_5   = create_float(0.5)
    float_0_25  = create_float(0.25)

    # Initialise return code.
    code = [""]

    # Create zero array for basisvalues.
    # Get number of members of the expansion set.
    num_mem = data["num_expansion_members"]
    code += [f_comment("Array of basisvalues.")]
    code += [f_decl(f_double, f_component(f_basisvalue, num_mem), f_tensor([0.0]*num_mem))]

    # Declare helper variables, will be removed if not used.
    code += ["", f_comment("Declare helper variables.")]
    code += [f_decl(f_uint, idx0, int_0)]
    code += [f_decl(f_uint, idx1, int_0)]
    code += [f_decl(f_uint, idx2, int_0)]
    code += [f_decl(f_double, str(an), f_float(0))]
    code += [f_decl(f_double, str(bn), f_float(0))]
    code += [f_decl(f_double, str(cn), f_float(0))]

    # Get the element cell domain.
    # FIXME: KBO: Change this when supporting R^2 in R^3 elements.
    element_cell_domain = data["cell_domain"]

    # 1D
    if (element_cell_domain == "interval"):
        # FIAT_NEW.expansions.LineExpansionSet.
        # FIAT_NEW code
        # psitilde_as = jacobi.eval_jacobi_batch(0,0,n,ref_pts)
        # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs)
        # The initial value basisvalue 0 is always 1.0
        # FIAT_NEW code
        # for ii in range(result.shape[1]):
        #    result[0,ii] = 1.0 + xs[ii,0] - xs[ii,0]
        code += ["", f_comment("Compute basisvalues.")]
        code += [f_assign(f_component(f_basisvalue, 0), f_float(1.0))]

        # Only continue if the embedded degree is larger than zero.
        if embedded_degree > 0:

            # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs).
            # result[1,:] = 0.5 * ( a - b + ( a + b + 2.0 ) * xsnew )
            # The initial value basisvalue 1 is always x
            code += [f_assign(f_component(f_basisvalue, 1), f_x)]

            # Only active is embedded_degree > 1.
            if embedded_degree > 1:
                # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs).
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
                # Declare helper variables (Note the a2 is always zero and therefore left out).
                code += [f_decl(f_double, str(f1), f_float(0))]
                code += [f_decl(f_double, str(f2), f_float(0))]
                code += [f_decl(f_double, str(f3), f_float(0))]
                lines = []
                loop_vars = [(str(symbol_p), 2, int_n1)]
                # Create names.
                basis_k   = create_symbol(f_component(f_basisvalue, str(symbol_p)), CONST)
                basis_km1 = create_symbol(f_component(f_basisvalue, f_sub([str(symbol_p), int_1])), CONST)
                basis_km2 = create_symbol(f_component(f_basisvalue, f_sub([str(symbol_p), int_2])), CONST)
                # Compute helper variables.
                a1 = create_product([float_2, symbol_p, symbol_p, float_2*symbol_p - float_2])
                a3 = create_fraction(create_product([float_2*symbol_p,\
                                                     float_2*symbol_p - float_2,
                                                     float_2*symbol_p - float_1]), a1)
                a4 = create_fraction(create_product([float_4*symbol_p, symbol_p - float_1, symbol_p - float_1]), a1)
                lines.append(f_assign(str(f1), a1 ))
                lines.append(f_assign(str(f2), a3 ))
                lines.append(f_assign(str(f3), a4 ))
                # Compute value.
                lines.append(f_assign(str(basis_k), create_product([f2, symbol_x*basis_km1]) - f3*basis_km2))
                # Create loop (block of lines).
                code += f_loop(lines, loop_vars)

        # Scale values.
        # FIAT_NEW.expansions.LineExpansionSet.
        # FIAT_NEW code
        # results = numpy.zeros( ( n+1 , len(pts) ) , type( pts[0][0] ) )
        # for k in range( n + 1 ):
        #    results[k,:] = psitilde_as[k,:] * math.sqrt( k + 0.5 )
        lines = []
        loop_vars = [(str(symbol_p), 0, int_n1)]
        # Create names.
        basis_k = create_symbol(f_component(f_basisvalue, str(symbol_p)), CONST)
        # Compute value.
        fac1 = create_symbol( f_sqrt(str(symbol_p + float_0_5)), CONST )
        lines += [format["imul"](str(basis_k), str(fac1))]
        # Create loop (block of lines).
        code += f_loop(lines, loop_vars)
    # 2D
    elif (element_cell_domain == "triangle"):
        # FIAT_NEW.expansions.TriangleExpansionSet.

        # Compute helper factors
        # FIAT_NEW code
        # f1 = (1.0+2*x+y)/2.0
        # f2 = (1.0 - y) / 2.0
        # f3 = f2**2
        fac1 = create_fraction(float_1 + float_2*symbol_x + symbol_y, float_2)
        fac2 = create_fraction(float_1 - symbol_y, float_2)
        code += [f_decl(f_double, str(f1), fac1)]
        code += [f_decl(f_double, str(f2), fac2)]
        code += [f_decl(f_double, str(f3), f2*f2)]

        code += ["", f_comment("Compute basisvalues.")]
        # The initial value basisvalue 0 is always 1.0.
        # FIAT_NEW code
        # for ii in range( results.shape[1] ):
        #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
        code += [f_assign(f_component(f_basisvalue, 0), f_float(1.0))]

        # Only continue if the embedded degree is larger than zero.
        if embedded_degree > 0:

            # The initial value of basisfunction 1 is equal to f1.
            # FIAT_NEW code
            # results[idx(1,0),:] = f1
            code += [f_assign(f_component(f_basisvalue, 1), str(f1))]

            # NOTE: KBO: The order of the loops is VERY IMPORTANT!!
            # Only active is embedded_degree > 1.
            if embedded_degree > 1:
                # FIAT_NEW code (loop 1 in FIAT)
                # for p in range(1,n):
                #    a = (2.0*p+1)/(1.0+p)
                #    b = p / (p+1.0)
                #    results[idx(p+1,0)] = a * f1 * results[idx(p,0),:] \
                #        - p/(1.0+p) * f3 *results[idx(p-1,0),:]
                # FIXME: KBO: Is there an error in FIAT? why is b not used?
                lines = []
                loop_vars = [(str(symbol_p), 1, embedded_degree)]
                # Compute indices
                lines.append(f_assign(idx0, _idx2D(f_group(f_add([str(symbol_p), int_1])), int_0)))
                lines.append(f_assign(idx1, _idx2D(symbol_p, int_0)))
                lines.append(f_assign(idx2, _idx2D(f_group(f_sub([str(symbol_p), int_1])), int_0)))
                # Compute single helper variable an.
                lines.append(f_assign(str(an), create_fraction(float_2*symbol_p + float_1, symbol_p + float_1)))
                # Compute value.
                fac0 = create_product([an, f1, basis_idx1])
                fac1 = create_product([create_fraction(symbol_p, float_1 + symbol_p), f3, basis_idx2])
                lines.append(f_assign(str(basis_idx0), fac0 - fac1))
                # Create loop (block of lines).
                code += f_loop(lines, loop_vars)

            # FIAT_NEW code (loop 2 in FIAT).
            # for p in range(n):
            #    results[idx(p,1),:] = 0.5 * (1+2.0*p+(3.0+2.0*p)*y) \
            #        * results[idx(p,0)]
            lines = []
            loop_vars = [(str(symbol_p), 0, embedded_degree)]
            # Compute indices.
            lines.append(f_assign(idx0, _idx2D(symbol_p, int_1)))
            lines.append(f_assign(idx1, _idx2D(symbol_p, int_0)))
            # Compute value.
            fac0 = create_product([float_3 + float_2*symbol_p, symbol_y])
            fac1 = create_product([float_0_5, float_1 + float_2*symbol_p + fac0, basis_idx1])
            lines.append(f_assign(str(basis_idx0), fac1.expand().reduce_ops()))
            # Create loop (block of lines).
            code += f_loop(lines, loop_vars)

            # Only active is embedded_degree > 1.
            if embedded_degree > 1:
                # FIAT_NEW code (loop 3 in FIAT).
                # for p in range(n-1):
                #    for q in range(1,n-p):
                #        (a1,a2,a3) = jrc(2*p+1,0,q)
                #        results[idx(p,q+1),:] \
                #            = ( a1 * y + a2 ) * results[idx(p,q)] \
                #            - a3 * results[idx(p,q-1)]
                lines = []
                loop_vars = [(str(symbol_p), 0, int_nm1),\
                             (str(symbol_q), 1, f_sub([int_n, str(symbol_p)]))]
                # Compute indices.
                lines.append(f_assign(idx0, _idx2D(symbol_p, f_add([str(symbol_q), int_1]) )))
                lines.append(f_assign(idx1, _idx2D(symbol_p, symbol_q)))
                lines.append(f_assign(idx2, _idx2D(symbol_p, f_sub([str(symbol_q), int_1]))))
                # Comute all helper variables.
                jrc = _jrc(float_2*symbol_p + float_1, float_0, symbol_q)
                lines.append(f_assign(str(an), jrc[0]))
                lines.append(f_assign(str(bn), jrc[1]))
                lines.append(f_assign(str(cn), jrc[2]))
                # Compute value.
                fac0 = create_product([an * symbol_y + bn, basis_idx1])
                fac1 = cn * basis_idx2
                lines.append(f_assign(str(basis_idx0), fac0 - fac1))
                # Create loop (block of lines).
                code += f_loop(lines, loop_vars)

            # FIAT_NEW code (loop 4 in FIAT).
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        results[idx(p,q),:] *= math.sqrt((p+0.5)*(p+q+1.0))
            lines = []
            loop_vars = [(str(symbol_p), 0, int_n1), \
                         (str(symbol_q), 0, f_sub([int_n1, str(symbol_p)]))]
            # Compute indices.
            lines.append(f_assign(idx0, _idx2D(symbol_p, symbol_q)))
            # Compute value.
            fac0 = create_product([symbol_p + float_0_5, symbol_p + symbol_q + float_1])
            fac2 = create_symbol( f_sqrt(str(fac0)), CONST )
            lines += [f_imul(str(basis_idx0), fac2)]
            # Create loop (block of lines).
            code += f_loop(lines, loop_vars)

    # 3D
    elif (element_cell_domain == "tetrahedron"):
        # FIAT_NEW.expansions.TetrahedronExpansionSet.

        # Compute helper factors.
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
        code += [f_decl(f_double, str(f1), fac1)]
        code += [f_decl(f_double, str(f2), fac2)]
        code += [f_decl(f_double, str(f3), fac3)]
        code += [f_decl(f_double, str(f4), fac4)]
        code += [f_decl(f_double, str(f5), f4*f4)]

        code += ["", f_comment("Compute basisvalues.")]
        # The initial value basisvalue 0 is always 1.0.
        # FIAT_NEW code
        # for ii in range( results.shape[1] ):
        #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
        code += [f_assign(f_component(f_basisvalue, 0), f_float(1.0))]

        # Only continue if the embedded degree is larger than zero.
        if embedded_degree > 0:

            # The initial value of basisfunction 1 is equal to f1.
            # FIAT_NEW code
            # results[idx(1,0),:] = f1
            code += [f_assign(f_component(f_basisvalue, 1), str(f1))]

            # NOTE: KBO: The order of the loops is VERY IMPORTANT!!
            # Only active is embedded_degree > 1
            if embedded_degree > 1:

                # FIAT_NEW code (loop 1 in FIAT).
                # for p in range(1,n):
                #    a1 = ( 2.0 * p + 1.0 ) / ( p + 1.0 )
                #    a2 = p / (p + 1.0)
                #    results[idx(p+1,0,0)] = a1 * factor1 * results[idx(p,0,0)] \
                #        -a2 * factor2 * results[ idx(p-1,0,0) ]
                lines = []
                loop_vars = [(str(symbol_p), 1, embedded_degree)]
                # Compute indices.
                lines.append(f_assign(idx0, _idx3D(f_group(f_add([str(symbol_p), int_1])), int_0, int_0)))
                lines.append(f_assign(idx1, _idx3D(str(symbol_p)               , int_0, int_0)))
                lines.append(f_assign(idx2, _idx3D(f_group(f_sub([str(symbol_p), int_1])), int_0, int_0)))
                # Compute value.
                fac1 = create_fraction(float_2*symbol_p + float_1, symbol_p + float_1)
                fac2 = create_fraction(symbol_p, symbol_p + float_1)
                fac3 = create_product([fac1, f1, basis_idx1]) - create_product([fac2, f2, basis_idx2])
                lines.append(f_assign(str(basis_idx0), fac3))
                # Create loop (block of lines).
                code += f_loop(lines, loop_vars)

            # FIAT_NEW code (loop 2 in FIAT).
            # q = 1
            # for p in range(0,n):
            #    results[idx(p,1,0)] = results[idx(p,0,0)] \
            #        * ( p * (1.0 + y) + ( 2.0 + 3.0 * y + z ) / 2 )
            lines = []
            loop_vars = [(str(symbol_p), 0, embedded_degree)]
            # Compute indices.
            lines.append(f_assign(idx0, _idx3D(symbol_p, int_1, int_0)))
            lines.append(f_assign(idx1, _idx3D(symbol_p, int_0, int_0)))
            # Compute value.
            fac1 = create_fraction(float_2 + float_3*symbol_y + symbol_z, float_2)
            fac2 = create_product([symbol_p, float_1 + symbol_y])
            fac3 = create_product([basis_idx1, fac2 + fac1])
            lines.append(f_assign(str(basis_idx0), fac3))
            # Create loop (block of lines).
            code += f_loop(lines, loop_vars)

            # Only active is embedded_degree > 1.
            if embedded_degree > 1:
                # FIAT_NEW code (loop 3 in FIAT).
                # for p in range(0,n-1):
                #    for q in range(1,n-p):
                #        (aq,bq,cq) = jrc(2*p+1,0,q)
                #        qmcoeff = aq * factor3 + bq * factor4
                #        qm1coeff = cq * factor5
                #        results[idx(p,q+1,0)] = qmcoeff * results[idx(p,q,0)] \
                #            - qm1coeff * results[idx(p,q-1,0)]
                lines = []
                loop_vars = [(str(symbol_p), 0, int_nm1),\
                             (str(symbol_q), 1, f_sub([int_n, str(symbol_p)]))]
                # Compute indices.
                lines.append(f_assign(idx0, _idx3D(symbol_p, f_group(f_add([str(symbol_q), int_1])), int_0)))
                lines.append(f_assign(idx1, _idx3D(symbol_p, symbol_q                     , int_0)))
                lines.append(f_assign(idx2, _idx3D(symbol_p, f_group(f_sub([str(symbol_q), int_1])), int_0)))
                # Comute all helper variables.
                jrc = _jrc(float_2*symbol_p + float_1, float_0, symbol_q)
                lines.append(f_assign(str(an), jrc[0]))
                lines.append(f_assign(str(bn), jrc[1]))
                lines.append(f_assign(str(cn), jrc[2]))
                # Compute value.
                fac1 = create_product([an*f3 + bn*f4, basis_idx1]) - cn*f5*basis_idx2
                lines.append(f_assign(str(basis_idx0), fac1))
                # Create loop (block of lines).
                code += f_loop(lines, loop_vars)

            # FIAT_NEW code (loop 4 in FIAT).
            # now handle r=1
            # for p in range(n):
            #    for q in range(n-p):
            #        results[idx(p,q,1)] = results[idx(p,q,0)] \
            #            * ( 1.0 + p + q + ( 2.0 + q + p ) * z )
            lines = []
            loop_vars = [(str(symbol_p), 0, int_n),\
                         (str(symbol_q), 0, f_sub([int_n, str(symbol_p)]))]
            # Compute indices.
            lines.append(f_assign(idx0, _idx3D(symbol_p, symbol_q, int_1)))
            lines.append(f_assign(idx1, _idx3D(symbol_p, symbol_q, int_0)))
            # Compute value.
            fac1 = create_product([float_2 + symbol_p + symbol_q, symbol_z])
            fac2 = create_product([basis_idx1, float_1 + symbol_p + symbol_q + fac1])
            lines.append(f_assign(str(basis_idx0), fac2))
            # Create loop (block of lines).
            code += f_loop(lines, loop_vars)

            # Only active is embedded_degree > 1.
            if embedded_degree > 1:
                # FIAT_NEW code (loop 5 in FIAT).
                # general r by recurrence
                # for p in range(n-1):
                #     for q in range(0,n-p-1):
                #         for r in range(1,n-p-q):
                #             ar,br,cr = jrc(2*p+2*q+2,0,r)
                #             results[idx(p,q,r+1)] = \
                #                         (ar * z + br) * results[idx(p,q,r) ] \
                #                         - cr * results[idx(p,q,r-1) ]
                lines = []
                loop_vars = [(str(symbol_p), 0, int_nm1),\
                             (str(symbol_q), 0, f_sub([int_nm1, str(symbol_p)])),\
                             (str(symbol_r), 1, f_sub([int_n, str(symbol_p), str(symbol_q)]))]
                # Compute indices.
                lines.append(f_assign(idx0, _idx3D(symbol_p, symbol_q, f_add([str(symbol_r), int_1]))))
                lines.append(f_assign(idx1, _idx3D(symbol_p, symbol_q, symbol_r)))
                lines.append(f_assign(idx2, _idx3D(symbol_p, symbol_q, f_sub([str(symbol_r), int_1]))))
                # Comute all helper variables.
                jrc = _jrc(float_2*symbol_p + float_2*symbol_q + float_2, float_0, symbol_r)
                lines.append(f_assign(str(an), jrc[0]))
                lines.append(f_assign(str(bn), jrc[1]))
                lines.append(f_assign(str(cn), jrc[2]))
                # Compute value.
                fac1 = create_product([an*symbol_z + bn, basis_idx1]) - cn*basis_idx2
                lines.append(f_assign(str(basis_idx0), fac1))
                # Create loop (block of lines).
                code += f_loop(lines, loop_vars)

            # FIAT_NEW code (loop 6 in FIAT).
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        for r in range(n-p-q+1):
            #            results[idx(p,q,r)] *= math.sqrt((p+0.5)*(p+q+1.0)*(p+q+r+1.5))
            lines = []
            loop_vars = [(str(symbol_p), 0, int_n1),\
                         (str(symbol_q), 0, f_sub([int_n1, str(symbol_p)])),\
                         (str(symbol_r), 0, f_sub([int_n1, str(symbol_p), str(symbol_q)]))]
            # Compute indices.
            lines.append(f_assign(idx0, _idx3D(symbol_p, symbol_q, symbol_r)))
            # Compute value
            fac0 = create_product([symbol_p + float_0_5,\
                                   symbol_p + symbol_q + float_1,\
                                   symbol_p + symbol_q + symbol_r + float_1_5])
            fac2 = create_symbol( f_sqrt(str(fac0)), CONST )
            lines += [format["imul"](str(basis_idx0), fac2)]
            # Create loop (block of lines).
            code += f_loop(lines, loop_vars)
    else:
        error("Cannot compute basis values for shape: %d" % elemet_cell_domain)

    return code + [""]

