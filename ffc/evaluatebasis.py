"""Code generation for evaluation of finite element basis values. This
module generates code which is more or less a C++ representation of
the code found in FIAT."""

# Copyright (C) 2007-2010 Kristian B. Oelgaard
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2007-04-04
# Last changed: 2013-01-10
#
# Modified by Marie E. Rognes 2011
# Modified by Anders Logg 2013
#
# MER: The original module generated code that was more or less a C++
# representation of the code found in FIAT. I've modified this (for 2
# and 3D) to generate code that does the same as FIAT, but with loops
# unrolled etc, thus removing unnecessary computations at runtime.
# There might be some clean-ups required, specially after this.

# Python modules
import math
import numpy

# FFC modules
from ffc.log import error
from ffc.cpp import remove_unused, indent, format
from ffc.quadrature.symbolics import create_float, create_float, create_symbol,\
                                     create_product, create_sum, create_fraction, CONST

def _evaluate_basis_all(data):
    """Like evaluate_basis, but return the values of all basis functions (dofs)."""

    if isinstance(data, str):
        return format["exception"]("evaluate_basis_all: %s" % data)

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
    physical_value_size = data["physical_value_size"]
    space_dimension = data["space_dimension"]

    # Special case where space dimension is one (constant elements).
    if space_dimension == 1:
        code += [f_comment("Element is constant, calling evaluate_basis.")]
        code += [f_basis(format["int"](0), f_values)]
        return "\n".join(code)

    # Declare helper value to hold single dof values.
    code += [f_comment("Helper variable to hold values of a single dof.")]
    if physical_value_size == 1:
        code += [f_decl(f_double, f_dof_vals, f_float(0.0))]
    else:
        code += [f_decl(f_double,
                        f_component(f_dof_vals, physical_value_size),
                        f_tensor([0.0]*physical_value_size))]

    # Create loop over dofs that calls evaluate_basis for a single dof and
    # inserts the values into the global array.
    code += ["", f_comment("Loop dofs and call evaluate_basis")]
    lines_r = []
    loop_vars_r = [(f_r, 0, space_dimension)]

    if physical_value_size == 1:
        lines_r += [f_basis(f_r, f_ref_var(f_dof_vals))]
    else:
        lines_r += [f_basis(f_r, f_dof_vals)]

    if physical_value_size ==  1:
        lines_r += [f_assign(f_component(f_values, f_r), f_dof_vals)]
    else:
        index = format["matrix index"](f_r, f_s, physical_value_size)
        lines_s = [f_assign(f_component(f_values, index), f_component(f_dof_vals, f_s))]
        lines_r += f_loop(lines_s, [(f_s, 0, physical_value_size)])

    code += f_loop(lines_r, loop_vars_r)

    # Generate code (no need to remove unused).
    return "\n".join(code)

# From FIAT_NEW.polynomial_set.tabulate()
def _evaluate_basis(data):
    """Generate run time code to evaluate an element basisfunction at an
    arbitrary point. The value(s) of the basisfunction is/are
    computed as in FIAT as the dot product of the coefficients (computed at compile time)
    and basisvalues which are dependent on the coordinate and thus have to be computed at
    run time.

    The function should work for all elements supported by FIAT, but it remains
    untested for tensor valued elements."""

    if isinstance(data, str):
        return format["exception"]("evaluate_basis: %s" % data)

    # Prefetch formats.
    f_assign    = format["assign"]
    f_comment   = format["comment"]
    f_values    = format["argument values"]
    f_float     = format["floating point"]
    f_component = format["component"]

    # Initialise return code.
    code = []

    # Get the element cell name and geometric dimension.
    element_cellname = data["cellname"]
    gdim = data["geometric_dimension"]
    tdim = data["topological_dimension"]

    # Get code snippets for Jacobian, Inverse of Jacobian and mapping of
    # coordinates from physical element to the FIAT reference element.
    code += [format["compute_jacobian"](tdim, gdim)]
    code += [format["compute_jacobian_inverse"](tdim, gdim)]
    if data["needs_oriented"]:
        code += [format["orientation"](tdim, gdim)]
    code += ["", format["fiat coordinate map"](element_cellname, gdim)]

    # Get value shape and reset values. This should also work for TensorElement,
    # scalar are empty tuples, therefore (1,) in which case value_shape = 1.
    reference_value_size = data["reference_value_size"]
    code += ["", f_comment("Reset values")]
    if reference_value_size == 1:
        # Reset values as a pointer.
        code += [f_assign(format["dereference pointer"](f_values), f_float(0.0))]
    else:
        # Reset all values.
        code += [f_assign(f_component(f_values, i), f_float(0.0)) for i in range(reference_value_size)]

    # Create code for all basis values (dofs).
    dof_cases = []
    for dof in data["dof_data"]:
        dof_cases.append(_generate_dof_code(data, dof))
    code += [format["switch"](format["argument basis num"], dof_cases)]

    # Remove unused variables (from transformations and mappings) in code.
    code = remove_unused("\n".join(code))
    #code = "\n".join(code)
    return code

def _generate_dof_code(data, dof_data):
    """Generate code for a single basis element as the dot product of
    coefficients and basisvalues. Then apply transformation if applicable."""

    # Generate basisvalues.
    code = _compute_basisvalues(data, dof_data)

    # Tabulate coefficients.
    code += _tabulate_coefficients(dof_data)

    # Compute the value of the basisfunction as the dot product of the
    # coefficients and basisvalues and apply transformation.
    code += _compute_values(data, dof_data)

    return remove_unused("\n".join(code))

def _tabulate_coefficients(dof_data):
    """This function tabulates the element coefficients that are
    generated by FIAT at compile time."""

    # Prefetch formats to speed up code generation.
    f_comment       = format["comment"]
    f_table         = format["static const float declaration"]
    f_coefficients  = format["coefficients"]
    f_component     = format["component"]
    f_decl          = format["declaration"]
    f_tensor        = format["tabulate tensor"]
    f_new_line      = format["new line"]

    # Get coefficients from basis functions, computed by FIAT at compile time.
    coefficients = dof_data["coeffs"]

    # Initialise return code.
    code = [f_comment("Table(s) of coefficients")]

    # Get number of members of the expansion set.
    num_mem = dof_data["num_expansion_members"]

    # Generate tables for each component.
    for i, coeffs in enumerate(coefficients):

        # Varable name for coefficients.
        name = f_component(f_coefficients(i), num_mem)

        # Generate array of values.
        code += [f_decl(f_table, name, f_new_line + f_tensor(coeffs))] + [""]
    return code

def _compute_values(data, dof_data):
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
#    f_dof           = format["local dof"]
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

    tdim = data["topological_dimension"]
    gdim = data["geometric_dimension"]

    # Initialise return code.
    code = [f_comment("Compute value(s)")]

    # Get dof data.
    num_components = dof_data["num_components"]
    offset = dof_data["offset"]

    lines = []
    if data["reference_value_size"] != 1:
        # Loop number of components.
        for i in range(num_components):
            # Generate name and value to create matrix vector multiply.
            name = f_component(f_values, i + offset)
            value = f_mul([f_component(f_coefficients(i), f_r),\
                    f_component(f_basisvalues, f_r)])
            lines += [f_iadd(name, value)]
    else:
        # Generate name and value to create matrix vector multiply.
        name = f_deref_pointer(f_values)
        value = f_mul([f_component(f_coefficients(0), f_r),\
                f_component(f_basisvalues, f_r)])
        lines = [f_iadd(name, value)]

    # Get number of members of the expansion set.
    num_mem = dof_data["num_expansion_members"]
    loop_vars = [(f_r, 0, num_mem)]
    code += f_loop(lines, loop_vars)

    # Apply transformation if applicable.
    mapping = dof_data["mapping"]
    if mapping == "affine":
        pass
    elif mapping == "contravariant piola":
        code += ["", f_comment("Using contravariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))\
                  for i in range(num_components)]
        # Create names for inner product.
        basis_col = [f_tmp_ref(j) for j in range(tdim)]
        for i in range(gdim):
            # Create Jacobian.
            jacobian_row = [f_trans("J", i, j, gdim, tdim, None) for j in range(tdim)]

            # Create inner product and multiply by inverse of Jacobian.
            inner = f_group(f_inner(jacobian_row, basis_col))
            value = f_mul([f_inv(f_detJ(None)), inner])
            name = f_component(f_values, i + offset)
            code += [f_assign(name, value)]
    elif mapping == "covariant piola":
        code += ["", f_comment("Using covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))\
                  for i in range(num_components)]
        # Create names for inner product.
        tdim = data["topological_dimension"]
        gdim = data["geometric_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(tdim)]
        for i in range(gdim):
            # Create inverse of Jacobian.
            inv_jacobian_column = [f_trans("JINV", j, i, tdim, gdim, None) for j in range(tdim)]

            # Create inner product of basis values and inverse of Jacobian.
            value = f_group(f_inner(inv_jacobian_column, basis_col))
            name = f_component(f_values, i + offset)
            code += [f_assign(name, value)]
    else:
        error("Unknown mapping: %s" % mapping)

    return code

def _compute_basisvalues(data, dof_data):
    """From FIAT_NEW.expansions."""

    UNROLL = True

    # Prefetch formats to speed up code generation.
    f_comment     = format["comment"]
    f_add         = format["add"]
    f_mul         = format["mul"]
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
    embedded_degree = dof_data["embedded_degree"]

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
    num_mem = dof_data["num_expansion_members"]
    code += [f_comment("Array of basisvalues")]
    code += [f_decl(f_double, f_component(f_basisvalue, num_mem), f_tensor([0.0]*num_mem))]

    # Declare helper variables, will be removed if not used.
    code += ["", f_comment("Declare helper variables")]
    code += [f_decl(f_uint, idx0, int_0)]
    code += [f_decl(f_uint, idx1, int_0)]
    code += [f_decl(f_uint, idx2, int_0)]
    code += [f_decl(f_double, str(an), f_float(0))]
    code += [f_decl(f_double, str(bn), f_float(0))]
    code += [f_decl(f_double, str(cn), f_float(0))]

    # Get the element cell name
    element_cellname = data["cellname"]

    def _jrc(a, b, n):
        an = float( ( 2*n+1+a+b)*(2*n+2+a+b))/ float( 2*(n+1)*(n+1+a+b))
        bn = float( (a*a-b*b) * (2*n+1+a+b))/ float( 2*(n+1)*(2*n+a+b)*(n+1+a+b) )
        cn = float( (n+a)*(n+b)*(2*n+2+a+b))/ float( (n+1)*(n+1+a+b)*(2*n+a+b) )
        return (an,bn,cn)

    # 1D
    if (element_cellname == "interval"):
        # FIAT_NEW.expansions.LineExpansionSet.
        # FIAT_NEW code
        # psitilde_as = jacobi.eval_jacobi_batch(0,0,n,ref_pts)
        # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs)
        # The initial value basisvalue 0 is always 1.0
        # FIAT_NEW code
        # for ii in range(result.shape[1]):
        #    result[0,ii] = 1.0 + xs[ii,0] - xs[ii,0]
        code += ["", f_comment("Compute basisvalues")]
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

                # The below implements the above (with a = b = apb = 0)
                for r in range(2, embedded_degree+1):

                    # Define helper variables
                    a1 = 2.0*r*r*(2.0*r - 2.0)
                    a3 = ((2.0*r - 2.0)*(2.0*r - 1.0 )*(2.0*r))/a1
                    a4 = (2.0*(r - 1.0)*(r - 1.0)*(2.0*r))/a1

                    assign_to = f_component(f_basisvalue, r)
                    assign_from = f_sub([f_mul([f_x, f_component(f_basisvalue, r-1), f_float(a3)]),
                                         f_mul([f_component(f_basisvalue, r-2), f_float(a4)])])
                    code += [f_assign(assign_to, assign_from)]

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
    elif (element_cellname == "triangle"):
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

        code += ["", f_comment("Compute basisvalues")]
        # The initial value basisvalue 0 is always 1.0.
        # FIAT_NEW code
        # for ii in range( results.shape[1] ):
        #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
        code += [f_assign(f_component(f_basisvalue, 0), f_float(1.0))]

        def _idx2d(p, q):
            return (p + q)*(p + q + 1)//2 + q

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

                for r in range(1, embedded_degree):
                    rr = _idx2d((r + 1), 0)
                    assign_to = f_component(f_basisvalue, rr)
                    ss = _idx2d(r, 0)
                    tt = _idx2d((r - 1), 0)
                    A = (2*r + 1.0)/(r + 1)
                    B = r/(1.0 + r)
                    v1 = f_mul([f_component(f_basisvalue, ss), f_float(A),
                                str(f1)])
                    v2 = f_mul([f_component(f_basisvalue, tt), f_float(B),
                                str(f3)])
                    assign_from = f_sub([v1, v2])
                    code += [f_assign(assign_to, assign_from)]

            # FIAT_NEW code (loop 2 in FIAT).
            # for p in range(n):
            #    results[idx(p,1),:] = 0.5 * (1+2.0*p+(3.0+2.0*p)*y) \
            #        * results[idx(p,0)]

            for r in range(0, embedded_degree):
                # (p+q)*(p+q+1)/2 + q
                rr = _idx2d(r, 1)
                assign_to = f_component(f_basisvalue, rr)
                ss = _idx2d(r, 0)
                A = 0.5*(1 + 2*r)
                B = 0.5*(3 + 2*r)
                C = f_add([f_float(A), f_mul([f_float(B), str(symbol_y)])])
                assign_from = f_mul([f_component(f_basisvalue, ss),
                                     f_group(C)])
                code += [f_assign(assign_to, assign_from)]

            # Only active is embedded_degree > 1.
            if embedded_degree > 1:
                # FIAT_NEW code (loop 3 in FIAT).
                # for p in range(n-1):
                #    for q in range(1,n-p):
                #        (a1,a2,a3) = jrc(2*p+1,0,q)
                #        results[idx(p,q+1),:] \
                #            = ( a1 * y + a2 ) * results[idx(p,q)] \
                #            - a3 * results[idx(p,q-1)]

                for r in range(0, embedded_degree - 1):
                    for s in range(1, embedded_degree - r):
                        rr = _idx2d(r, (s + 1))
                        ss = _idx2d(r, s)
                        tt = _idx2d(r, s - 1)
                        A, B, C = _jrc(2*r + 1, 0, s)
                        assign_to = f_component(f_basisvalue, rr)
                        assign_from = f_sub([f_mul([f_component(f_basisvalue, ss), f_group(f_add([f_float(B), f_mul([str(symbol_y), f_float(A)])]))]),
                                             f_mul([f_component(f_basisvalue, tt), f_float(C)])])
                        code += [f_assign(assign_to, assign_from)]

            # FIAT_NEW code (loop 4 in FIAT).
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        results[idx(p,q),:] *= math.sqrt((p+0.5)*(p+q+1.0))
            n1 = embedded_degree + 1
            for r in range(0, n1):
                for s in range(0, n1 - r):
                    rr = _idx2d(r, s)
                    A = (r + 0.5)*(r + s + 1)
                    assign_to = f_component(f_basisvalue, rr)
                    code += [f_imul(assign_to, f_sqrt(A))]

    # 3D
    elif (element_cellname == "tetrahedron"):

        # FIAT_NEW code (compute index function) TetrahedronExpansionSet.
        # def idx(p,q,r):
        #     return (p+q+r)*(p+q+r+1)*(p+q+r+2)/6 + (q+r)*(q+r+1)/2 + r
        def _idx3d(p, q, r):
            return (p+q+r)*(p+q+r+1)*(p+q+r+2)/6 + (q+r)*(q+r+1)/2 + r

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

        code += ["", f_comment("Compute basisvalues")]
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
                for r in range(1, embedded_degree):
                    rr = _idx3d((r + 1), 0, 0)
                    ss = _idx3d(r, 0, 0)
                    tt = _idx3d((r - 1), 0, 0)
                    A = (2*r + 1.0)/(r + 1)
                    B = r/(r + 1.0)
                    assign_to = f_component(f_basisvalue, rr)
                    assign_from = f_sub([f_mul([f_float(A), str(f1), f_component(f_basisvalue, ss)]), f_mul([f_float(B), str(f2), f_component(f_basisvalue, tt)])])
                    code += [f_assign(assign_to, assign_from)]

            # FIAT_NEW code (loop 2 in FIAT).
            # q = 1
            # for p in range(0,n):
            #    results[idx(p,1,0)] = results[idx(p,0,0)] \
            #        * ( p * (1.0 + y) + ( 2.0 + 3.0 * y + z ) / 2 )

            for r in range(0, embedded_degree):
                rr = _idx3d(r, 1, 0)
                ss = _idx3d(r, 0, 0)
                assign_to = f_component(f_basisvalue, rr)
                term0 = f_mul([f_float(0.5), f_group(f_add([f_float(2), f_mul([f_float(3), str(symbol_y)]), str(symbol_z)]))])
                if r == 0:
                    assign_from = f_mul([term0, f_component(f_basisvalue, ss)])
                else:
                    term1 = f_mul([f_float(r), f_group(f_add([f_float(1), str(symbol_y)]))])
                    assign_from = f_mul([f_group(f_add([term0, term1])), f_component(f_basisvalue, ss)])
                code += [f_assign(assign_to, assign_from)]

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

                for r in range(0, embedded_degree - 1):
                    for s in range(1, embedded_degree - r):
                        rr = _idx3d(r, (s + 1), 0)
                        ss = _idx3d(r, s, 0)
                        tt = _idx3d(r, s - 1, 0)
                        (A, B, C) = _jrc(2*r + 1, 0, s)
                        assign_to = f_component(f_basisvalue, rr)
                        term0 = f_mul([f_group(f_add([f_mul([f_float(A), str(f3)]), f_mul([f_float(B), str(f4)])])), f_component(f_basisvalue, ss)])
                        term1 = f_mul([f_float(C), str(f5), f_component(f_basisvalue, tt)])
                        assign_from = f_sub([term0, term1])
                        code += [f_assign(assign_to, assign_from)]

            # FIAT_NEW code (loop 4 in FIAT).
            # now handle r=1
            # for p in range(n):
            #    for q in range(n-p):
            #        results[idx(p,q,1)] = results[idx(p,q,0)] \
            #            * ( 1.0 + p + q + ( 2.0 + q + p ) * z )
            for r in range(0, embedded_degree):
                for s in range(0, embedded_degree - r):
                    rr = _idx3d(r, s, 1)
                    ss = _idx3d(r, s, 0)
                    assign_to = f_component(f_basisvalue, rr)
                    A = f_add([f_mul([f_float(2 + r + s), str(symbol_z)]), f_float(1 + r + s)])
                    assign_from = f_mul([f_group(A), f_component(f_basisvalue, ss)])
                    code += [f_assign(assign_to, assign_from)]

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
                for r in range(embedded_degree - 1):
                    for s in range(0, embedded_degree - r - 1):
                        for t in range(1, embedded_degree - r - s):
                            rr = _idx3d(r, s, ( t + 1))
                            ss = _idx3d(r, s, t)
                            tt = _idx3d(r, s, t - 1)

                            (A, B, C) = _jrc(2*r + 2*s + 2, 0, t)
                            assign_to = f_component(f_basisvalue, rr)
                            az_b = f_group(f_add([f_float(B), f_mul([f_float(A), str(symbol_z)])]))
                            assign_from = f_sub([f_mul([f_component(f_basisvalue, ss), az_b]), f_mul([f_float(C), f_component(f_basisvalue, tt)])])
                            code += [f_assign(assign_to, assign_from)]

            # FIAT_NEW code (loop 6 in FIAT).
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        for r in range(n-p-q+1):
            #            results[idx(p,q,r)] *= math.sqrt((p+0.5)*(p+q+1.0)*(p+q+r+1.5))
            for r in range(embedded_degree + 1):
                for s in range(embedded_degree - r + 1):
                    for t in range(embedded_degree - r - s + 1):
                        rr = _idx3d(r, s, t)
                        A = (r + 0.5)*(r + s + 1)*(r + s + t + 1.5)
                        assign_to = f_component(f_basisvalue, rr)
                        multiply_by = f_sqrt(A)
                        myline = f_imul(assign_to, multiply_by)
                        code += [myline]

    else:
        error("Cannot compute basis values for shape: %d" % elemet_cell_domain)

    return code + [""]
