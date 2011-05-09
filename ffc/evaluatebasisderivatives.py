"""Code generation for evaluation of derivatives of finite element basis values.
This module generates code which is more or less a C++ representation of the code
found in FIAT_NEW."""

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC.  If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2007-04-16
# Last changed: 2011-02-21

# Python modules
import math
import numpy

# FFC modules
from ffc.log import error, ffc_assert
from ffc.evaluatebasis import _compute_basisvalues, _tabulate_coefficients
from ffc.cpp import remove_unused, indent, format

def _evaluate_basis_derivatives_all(data):
    """Like evaluate_basis, but return the values of all basis functions (dofs)."""

    if isinstance(data, str):
        return format["exception"]("evaluate_basis_derivatives_all: %s" % data)

    # Prefetch formats.
    f_r, f_s      = format["free indices"][:2]
    f_assign      = format["assign"]
    f_loop        = format["generate loop"]
    f_array       = format["dynamic array"]
    f_dof_vals    = format["dof values"]
    f_comment     = format["comment"]
    f_derivs      = format["call basis_derivatives"]
    f_values      = format["argument values"]
    f_int         = format["int"]
    f_num_derivs  = format["num derivatives"]
    f_double      = format["float declaration"]
    f_component   = format["component"]
    f_mul         = format["mul"]
    f_float       = format["floating point"]
    f_index       = format["matrix index"]
    f_del_array   = format["delete dynamic array"]

    # Initialise return code
    code = []

    # FIXME: KBO: Figure out which return format to use, either:
    # [dN0[0]/dx, dN0[0]/dy, dN0[1]/dx, dN0[1]/dy, dN1[0]/dx, dN1[0]/dy, dN1[1]/dx, dN1[1]/dy, ...]
    # or
    # [dN0[0]/dx, dN1[0]/dx, ..., dN0[1]/dx, dN1[1]/dx, ..., dN0[0]/dy, dN1[0]/dy, ..., dN0[1]/dy, dN1[1]/dy, ...]
    # or
    # [dN0[0]/dx, dN0[1]/dx, ..., dN1[0]/dx, dN1[1]/dx, ..., dN0[0]/dy, dN0[1]/dy, ..., dN1[0]/dy, dN1[1]/dy, ...]
    # for vector (tensor elements), currently returning option 1.

    # FIXME: KBO: For now, just call evaluate_basis_derivatives and map values
    # accordingly, this will keep the amount of code at a minimum. If it turns
    # out that speed is an issue (overhead from calling evaluate_basis), we can
    # easily generate all the code.

    # Get total value shape and space dimension for entire element (possibly mixed).
    value_size = data["value_size"]
    space_dimension = data["space_dimension"]

    # Special case where space dimension is one (constant elements).
    if space_dimension == 1:
        code += [f_comment("Element is constant, calling evaluate_basis_derivatives.")]
        code += [f_derivs(f_int(0), f_values)]
        return "\n".join(code)

    # Compute number of derivatives.
    code += _compute_num_derivatives(data["topological_dimension"])

    # Declare helper value to hold single dof values and reset.
    code += ["", f_comment("Helper variable to hold values of a single dof.")]
    if (value_size == 1):
        num_vals = f_num_derivs
    else:
        num_vals = f_mul([f_int(value_size), f_num_derivs])
    code += [f_array(f_double, f_dof_vals, num_vals)]
    line  = [f_assign(f_component(f_dof_vals, f_r), f_float(0.0))]
    code += f_loop(line, [(f_r, 0, num_vals)])

    # Create loop over dofs that calls evaluate_basis_derivatives for a single dof and
    # inserts the values into the global array.
    code += ["", f_comment("Loop dofs and call evaluate_basis_derivatives.")]
    name  = f_component(f_values, f_index(f_r, f_s, num_vals))
    value = f_component(f_dof_vals, f_s)
    lines_s  = [f_assign(name, value)]
    loop_s   = [(f_s, 0, num_vals)]

    lines_r  = [f_derivs(f_r, f_dof_vals)]
    lines_r += f_loop(lines_s, loop_s)
    loop_r   = [(f_r, 0, space_dimension)]
    code    += f_loop(lines_r, loop_r)

    code += ["", f_comment("Delete pointer.")]
    code += [f_del_array(f_dof_vals)]

    # Generate bode (no need to remove unused).
    return "\n".join(code)

def _evaluate_basis_derivatives(data):
    """Evaluate the derivatives of an element basisfunction at a point. The values are
    computed as in FIAT as the matrix product of the coefficients (computed at compile time),
    basisvalues which are dependent on the coordinate and thus have to be computed at
    run time and combinations (depending on the order of derivative) of dmats
    tables which hold the derivatives of the expansion coefficients."""

    if isinstance(data, str):
        return format["exception"]("evaluate_basis_derivatives: %s" % data)

    # Initialise return code.
    code = []

    # Get the element cell domain, geometric and topological dimension.
    element_cell_domain = data["cell_domain"]
    geometric_dimension = data["geometric_dimension"]
    topological_dimension = data["topological_dimension"]

    # Get code snippets for Jacobian, Inverse of Jacobian and mapping of
    # coordinates from physical element to the FIAT reference element.
    # FIXME: KBO: Change this when supporting R^2 in R^3 elements.
    code += [format["jacobian and inverse"](geometric_dimension)]
    code += ["", format["fiat coordinate map"](element_cell_domain)]

    # Compute number of derivatives that has to be computed, and declare an array to hold
    # the values of the derivatives on the reference element.
    code += [""]
    code += _compute_num_derivatives(topological_dimension)

    # Generate all possible combinations of derivatives.
    code += _generate_combinations(topological_dimension)

    # Generate the transformation matrix.
    code += _generate_transform(element_cell_domain)

    # Reset all values.
    code += _reset_values(data)

    # Create code for all basis values (dofs).
    dof_cases = []
    for dof in data["dof_data"]:
        dof_cases.append(_generate_dof_code(data, dof))
    code += [format["switch"](format["argument basis num"], dof_cases)]
    code = remove_unused("\n".join(code))
#    code = "\n".join(code)
    return code

def _compute_num_derivatives(topological_dimension):
    "Computes the number of derivatives of order 'n' as: element.topological_dimension()^n."
    # Prefetch formats.
    f_int         = format["int"]
    f_num_derivs  = format["num derivatives"]

    # Use loop to compute power since using std::pow() result in an ambiguous call.
    code = [format["comment"]("Compute number of derivatives.")]
    code.append(format["declaration"](format["uint declaration"], f_num_derivs, f_int(1)))
    loop_vars = [(format["free indices"][0], 0, format["argument derivative order"])]
    lines = [format["imul"](f_num_derivs, f_int(topological_dimension))]
    code += format["generate loop"](lines, loop_vars)

    return code

def _generate_combinations(topological_dimension):
    "Generate all possible combinations of derivatives of order 'n'."

    # Use code from format.
    code = ["", format["combinations"]\
            % {"combinations": format["derivative combinations"],\
               "topological_dimension-1": topological_dimension-1,\
               "num_derivatives" : format["num derivatives"],\
               "n": format["argument derivative order"]}]
    return code

def _generate_transform(element_cell_domain):
    """Generate the transformation matrix, whic is used to transform derivatives from reference
    element back to the physical element."""

    # Generate code to construct the inverse of the Jacobian
    if (element_cell_domain in ["interval", "triangle", "tetrahedron"]):
        code = ["", format["transform snippet"][element_cell_domain]\
        % {"transform": format["transform matrix"],\
           "num_derivatives" : format["num derivatives"],\
           "n": format["argument derivative order"],\
           "combinations": format["derivative combinations"],\
           "K":format["transform Jinv"]}]
    else:
        error("Cannot generate transform for shape: %s" % element_cell_domain)

    return code

def _reset_values(data):
    "Reset all components of the 'values' array as it is a pointer to an array."

    # Prefetch formats.
    f_assign  = format["assign"]
    f_r       = format["free indices"][0]

    code = ["", format["comment"]("Reset values. Assuming that values is always an array.")]

    # Get value shape and reset values. This should also work for TensorElement,
    # scalar are empty tuples, therefore (1,) in which case value_shape = 1.
    value_size = data["value_size"]

    # Only multiply by value shape if different from 1.
    if value_size == 1:
        num_vals = format["num derivatives"]
    else:
        num_vals = format["mul"]([format["int"](value_size), format["num derivatives"]])
    name = format["component"](format["argument values"], f_r)
    loop_vars = [(f_r, 0, num_vals)]
    lines = [f_assign(name, format["floating point"](0))]
    code += format["generate loop"](lines, loop_vars)

    return code + [""]

def _generate_dof_code(data, dof_data):
    "Generate code for a basis."

    code = []

    # Compute basisvalues, from evaluatebasis.py.
    code += _compute_basisvalues(data, dof_data)

    # Tabulate coefficients.
    code += _tabulate_coefficients(dof_data)

    # Tabulate coefficients for derivatives.
    code += _tabulate_dmats(dof_data)

    # Compute the derivatives of the basisfunctions on the reference (FIAT) element,
    # as the dot product of the new coefficients and basisvalues.
    code += _compute_reference_derivatives(data, dof_data)

    # Transform derivatives to physical element by multiplication with the transformation matrix.
    code += _transform_derivatives(data, dof_data)

    # Delete pointers.
    code += _delete_pointers()

    code = remove_unused("\n".join(code))
#    code = "\n".join(code)

    return code

#def _mixed_elements(data_list):
#    "Generate code for each sub-element in the event of vector valued elements or mixed elements."

#    # Prefetch formats to speed up code generation.
#    f_dofmap_if = format["dofmap if"]
#    f_if        = format["if"]

#    sum_value_dim = 0
#    sum_space_dim = 0

#    # Init return code.
#    code = []

#    # Generate code for each element.
#    for data in data_list:

#        # Get value and space dimension (should be tensor ready).
#        value_dim = sum(data["value_shape"] or (1,))
#        space_dim = data["space_dimension"]

#        # Generate map from global to local dof.
#        element_code = [_map_dof(sum_space_dim)]

#        # Generate code for basis element.
#        element_code += _generate_element_code(data, sum_value_dim)

#        # Remove unused code for each sub element and indent code.
#        if_code = indent(remove_unused("\n".join(element_code)), 2)

#        # Create if statement and add to code.
#        code += [f_if(f_dofmap_if(sum_space_dim, sum_space_dim + space_dim - 1), if_code)]

#        # Increase sum of value dimension, and space dimension.
#        sum_value_dim += value_dim
#        sum_space_dim += space_dim

#    return code

def _tabulate_dmats(dof_data):
    "Tabulate the derivatives of the polynomial base"

    code = []

    # Prefetch formats to speed up code generation.
    f_table     = format["static const float declaration"]
    f_dmats     = format["dmats"]
    f_component = format["component"]
    f_decl      = format["declaration"]
    f_tensor    = format["tabulate tensor"]
    f_new_line  = format["new line"]

    # Get derivative matrices (coefficients) of basis functions, computed by FIAT at compile time.
    derivative_matrices = dof_data["dmats"]

    code += [format["comment"]("Tables of derivatives of the polynomial base (transpose).")]

    # Generate tables for each spatial direction.
    for i, dmat in enumerate(derivative_matrices):

        # Extract derivatives for current direction (take transpose, FIAT_NEW PolynomialSet.tabulate()).
        matrix = numpy.transpose(dmat)

        # Get shape and check dimension (This is probably not needed).
        shape = numpy.shape(matrix)
        ffc_assert(shape[0] == shape[1] == dof_data["num_expansion_members"], "Something is wrong with the shape of dmats.")

        # Declare varable name for coefficients.
        name = f_component(f_dmats(i), [shape[0], shape[1]])
        code += [f_decl(f_table, name, f_new_line + f_tensor(matrix)), ""]

    return code

def _reset_dmats(shape_dmats, indices):
    "Set values in dmats equal to the identity matrix."
    f_assign  = format["assign"]
    f_float   = format["floating point"]
    i,j = indices

    code = [format["comment"]("Resetting dmats values to compute next derivative.")]
    dmats_old = format["component"](format["dmats"](""), [i, j])
    lines = [f_assign(dmats_old, f_float(0.0))]
    lines += [format["if"](i + format["is equal"] + j,\
              f_assign(dmats_old, f_float(1.0)))]
    loop_vars = [(i, 0, shape_dmats[0]), (j, 0, shape_dmats[1])]
    code += format["generate loop"](lines, loop_vars)
    return code

def _update_dmats(shape_dmats, indices):
    "Update values in dmats_old with values in dmats and set values in dmats to zero."
    f_assign    = format["assign"]
    f_component = format["component"]
    i,j = indices

    code = [format["comment"]("Updating dmats_old with new values and resetting dmats.")]
    dmats = f_component(format["dmats"](""), [i, j])
    dmats_old = f_component(format["dmats old"], [i, j])
    lines = [f_assign(dmats_old, dmats), f_assign(dmats, format["floating point"](0.0))]
    loop_vars = [(i, 0, shape_dmats[0]), (j, 0, shape_dmats[1])]
    code += format["generate loop"](lines, loop_vars)
    return code

def _compute_dmats(num_dmats, shape_dmats, available_indices, deriv_index):
    "Compute values of dmats as a matrix product."
    f_comment = format["comment"]
    s, t, u = available_indices

    # Reset dmats_old
    code = _reset_dmats(shape_dmats, [t, u])
    code += ["", f_comment("Looping derivative order to generate dmats.")]

    # Set dmats matrix equal to dmats_old
    lines = _update_dmats(shape_dmats, [t, u])

    lines += ["", f_comment("Update dmats using an inner product.")]
    # Create dmats matrix by multiplication
    comb = format["component"](format["derivative combinations"], [deriv_index, s])
    for i in range(num_dmats):
        lines += _dmats_product(shape_dmats, comb, i, [t, u])

    loop_vars = [(s, 0, format["argument derivative order"])]
    code += format["generate loop"](lines, loop_vars)

    return code

def _dmats_product(shape_dmats, index, i, indices):
    "Create product to update dmats."
    f_loop      = format["generate loop"]
    f_component = format["component"]
    t, u = indices
    tu = t + u

    dmats = f_component(format["dmats"](""), [t, u])
    dmats_old = f_component(format["dmats old"], [tu, u])
    value = format["multiply"]([f_component(format["dmats"](i), [t, tu]), dmats_old])
    name = format["iadd"](dmats, value)
    lines = f_loop([name], [(tu, 0, shape_dmats[0])])
    loop_vars = [(t, 0, shape_dmats[0]), (u, 0, shape_dmats[1])]
    code = [format["if"](index + format["is equal"] + str(i),\
            "\n".join(f_loop(lines, loop_vars)))]

    return code

def _compute_reference_derivatives(data, dof_data):
    """Compute derivatives on the reference element by recursively multiply coefficients with
    the relevant derivatives of the polynomial base until the requested order of derivatives
    has been reached. After this take the dot product with the basisvalues."""

    # Prefetch formats to speed up code generation
    f_comment       = format["comment"]
    f_num_derivs    = format["num derivatives"]
    f_mul           = format["mul"]
    f_int           = format["int"]
    f_matrix_index  = format["matrix index"]
    f_coefficients  = format["coefficients"]
#    f_dof           = format["local dof"]
    f_basisvalues   = format["basisvalues"]
    f_const_double  = format["const float declaration"]
    f_group         = format["grouping"]
    f_transform     = format["transform"]
    f_double        = format["float declaration"]
    f_component     = format["component"]
    f_tmp           = format["tmp ref value"]
    f_dmats         = format["dmats"]
    f_dmats_old     = format["dmats old"]
    f_assign        = format["assign"]
    f_decl          = format["declaration"]
    f_iadd          = format["iadd"]
    f_add           = format["add"]
    f_tensor        = format["tabulate tensor"]
    f_new_line      = format["new line"]
    f_loop          = format["generate loop"]
    f_derivatives   = format["reference derivatives"]
    f_array         = format["dynamic array"]
    f_float         = format["floating point"]
    f_inv           = format["inverse"]
    f_detJ          = format["det(J)"]

    f_r, f_s, f_t, f_u = format["free indices"]

    # Get number of components.
    num_components = dof_data["num_components"]

    # Get shape of derivative matrix (they should all have the same shape) and
    # verify that it is a square matrix.
    shape_dmats = numpy.shape(dof_data["dmats"][0])
    ffc_assert(shape_dmats[0] == shape_dmats[1],\
               "Something is wrong with the dmats:\n%s" % str(dof_data["dmats"]))

    code = [f_comment("Compute reference derivatives.")]

    # Declare pointer to array that holds derivatives on the FIAT element
    code += [f_comment("Declare pointer to array of derivatives on FIAT element.")]
    # The size of the array of reference derivatives is equal to the number of derivatives
    # times the number of components of the basis element
    if (num_components == 1):
        num_vals = f_num_derivs
    else:
        num_vals = f_mul([f_int(num_components), f_num_derivs])
    code += [f_array(f_double, f_derivatives, num_vals)]
    # Reset values of reference derivatives.
    name = f_component(f_derivatives, f_r)
    lines = [f_assign(name, f_float(0))]
    code += f_loop(lines, [(f_r, 0, num_vals)])
    code += [""]

    # Declare matrix of dmats (which will hold the matrix product of all combinations)
    # and dmats_old which is needed in order to perform the matrix product.
    value = f_tensor(numpy.eye(shape_dmats[0]))
    code += [f_comment("Declare derivative matrix (of polynomial basis).")]
    name = f_component(f_dmats(""), [shape_dmats[0], shape_dmats[1]])
    code += [f_decl(f_double, name, f_new_line + value), ""]
    code += [f_comment("Declare (auxiliary) derivative matrix (of polynomial basis).")]
    name = f_component(f_dmats_old, [shape_dmats[0], shape_dmats[1]])
    code += [f_decl(f_double, name, f_new_line + value), ""]

    # Compute dmats as a recursive matrix product
    lines = _compute_dmats(len(dof_data["dmats"]), shape_dmats, [f_s, f_t, f_u], f_r)

    # Compute derivatives for all components
    lines_c = []
    for i in range(num_components):
        name = f_component(f_derivatives, f_matrix_index(i, f_r, f_num_derivs))
        coeffs = f_component(f_coefficients(i), f_s)
        dmats = f_component(f_dmats(""), [f_s, f_t])
        basis = f_component(f_basisvalues, f_t)
        lines_c.append(f_iadd(name, f_mul([coeffs, dmats, basis])))
    loop_vars_c = [(f_s, 0, shape_dmats[0]),(f_t, 0, shape_dmats[1])]
    lines += f_loop(lines_c, loop_vars_c)

    # Apply transformation if applicable.
    mapping = dof_data["mapping"]
    if mapping == "affine":
        pass
    elif mapping == "contravariant piola":
        lines += ["", f_comment\
                ("Using contravariant Piola transform to map values back to the physical element.")]
        # Get temporary values before mapping.
        lines += [f_const_double(f_tmp(i),\
                  f_component(f_derivatives, f_matrix_index(i, f_r, f_num_derivs))) for i in range(num_components)]

        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [f_tmp(j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create Jacobian.
            jacobian_row = [f_transform("J", i, j, None) for j in range(topological_dimension)]

            # Create inner product and multiply by inverse of Jacobian.
            inner = [f_mul([jacobian_row[j], basis_col[j]]) for j in range(topological_dimension)]
            sum_ = f_group(f_add(inner))
            value = f_mul([f_inv(f_detJ(None)), sum_])
            name = f_component(f_derivatives, f_matrix_index(i, f_r, f_num_derivs))
            lines += [f_assign(name, value)]
    elif mapping == "covariant piola":
        lines += ["", f_comment\
                ("Using covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        lines += [f_const_double(f_tmp(i),\
                  f_component(f_derivatives, f_matrix_index(i, f_r, f_num_derivs))) for i in range(num_components)]
        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [f_tmp(j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create inverse of Jacobian.
            inv_jacobian_column = [f_transform("JINV", j, i, None) for j in range(topological_dimension)]

            # Create inner product of basis values and inverse of Jacobian.
            inner = [f_mul([inv_jacobian_column[j], basis_col[j]]) for j in range(topological_dimension)]
            value = f_group(f_add(inner))
            name = f_component(f_derivatives, f_matrix_index(i, f_r, f_num_derivs))
            lines += [f_assign(name, value)]
    else:
        error("Unknown mapping: %s" % mapping)

    # Generate loop over number of derivatives.
    # Loop all derivatives and compute value of the derivative as:
    # deriv_on_ref[r] = coeff[dof][s]*dmat[s][t]*basis[t]
    code += [f_comment("Loop possible derivatives.")]
    loop_vars = [(f_r, 0, f_num_derivs)]
    code += f_loop(lines, loop_vars)

    return code + [""]

def _transform_derivatives(data, dof_data):
    """Transform derivatives back to the physical element by applying the
    transformation matrix."""

    # Prefetch formats to speed up code generation.
    f_loop        = format["generate loop"]
    f_num_derivs  = format["num derivatives"]
    f_derivatives = format["reference derivatives"]
    f_values      = format["argument values"]
    f_mul         = format["mul"]
    f_iadd        = format["iadd"]
    f_component   = format["component"]
    f_transform   = format["transform matrix"]
    f_r, f_s      = format["free indices"][:2]
    f_index       = format["matrix index"]

    # Get number of components and offset.
    num_components = dof_data["num_components"]
    offset = dof_data["offset"]

    code = [format["comment"]("Transform derivatives back to physical element")]

    lines = []
    for i in range(num_components):
        access_name = f_index(offset + i, f_r, f_num_derivs)
        name = f_component(f_values, access_name)
        access_val = f_index(i, f_s, f_num_derivs)
        value = f_mul([f_component(f_transform, [f_r, f_s]), f_component(f_derivatives, access_val)])
        lines += [f_iadd(name, value)]

    loop_vars = [(f_r, 0, f_num_derivs), (f_s, 0, f_num_derivs)]
    code += f_loop(lines, loop_vars)
    return code

def _delete_pointers():
    "Delete the pointers to arrays."

    f_del_array = format["delete dynamic array"]
    code = []

    code += ["", format["comment"]("Delete pointer to array of derivatives on FIAT element")]
    code += [f_del_array(format["reference derivatives"]), ""]

    code += [format["comment"]("Delete pointer to array of combinations of derivatives and transform")]
    code += [f_del_array(format["derivative combinations"], format["num derivatives"])]
    code += [f_del_array(format["transform matrix"], format["num derivatives"])]

    return code

