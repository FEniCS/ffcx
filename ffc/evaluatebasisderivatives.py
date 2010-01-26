"""Code generation for evaluation of derivatives of finite element basis values.
This module generates code which is more or less a C++ representation of the code
found in FIAT_NEW."""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2007-04-16"
__copyright__ = "Copyright (C) 2007-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-25

# Python modules
import math
import numpy

# FFC modules
from ffc.log import error, ffc_assert
from ffc.evaluatebasis import _map_dof
from ffc.evaluatebasis import _compute_basisvalues
from ffc.evaluatebasis import _tabulate_coefficients
from ffc.cpp import tabulate_matrix, IndentControl, remove_unused, tabulate_vector
from ffc.quadrature.quadraturegenerator_utils import generate_loop
#from ffc.quadrature.symbolics import set_format

# Temporary import
#from cpp import format_old as format
from ffc.cpp import format

def _evaluate_basis_derivatives_all(data_list):
    """Like evaluate_basis, but return the values of all basis functions (dofs)."""

    format_r, format_s = format["free indices"][:2]
    format_assign = format["assign"]


    # Initialise objects
    Indent = IndentControl()
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
    value_shape = sum(sum(data["value_shape"] or (1,)) for data in data_list)
    space_dimension = sum(data["space_dimension"] for data in data_list)

    # Special case where space dimension is one (constant elements)
    if space_dimension == 1:
        code += [format["comment"]("Element is constant, calling evaluate_basis_derivatives.")]
        code += ["evaluate_basis_derivatives(0, n, %s, coordinates, c);" % format["argument values"]]
        return "\n".join(code)

    # Compute number of derivatives
    # Get the topological dimension.
    topological_dimension = data_list[0]["topological_dimension"]
    code += _compute_num_derivatives(topological_dimension, Indent, format)

    # Declare helper value to hold single dof values and reset
    code += ["", format["comment"]("Helper variable to hold values of a single dof.")]
    if (value_shape == 1):
        num_vals = format["num derivatives"]
    else:
        num_vals = format["multiply"]([format["floating point"](value_shape), format["num derivatives"]])
    code += [format_assign(Indent.indent(format["float declaration"] + format["pointer"] + "dof_values"),\
              format["component"](format["new"] + format["float declaration"], num_vals))]
    loop_vars = [(format_r, 0, num_vals)]
    line = [format_assign(format["component"]("dof_values", format_r), format["floating point"](0.0))]
    code += generate_loop(line, loop_vars, Indent, format)

    # Create loop over dofs that calls evaluate_basis_derivatives for a single dof and
    # inserts the values into the global array.
    code += ["", format["comment"]("Loop dofs and call evaluate_basis_derivatives.")]
    lines_r = []
    loop_vars_r = [(format_r, 0, space_dimension)]

    # FIXME: KBO: Move evaluate_basis_derivatives string to cpp.py
    lines_r += ["evaluate_basis_derivatives(%s, n, dof_values, coordinates, c);" % format_r]

    loop_vars_s = [(format_s, 0, num_vals)]
    index = format["add"]([format["multiply"]([format_r, num_vals]), format_s])
    name = format["component"](format["argument values"], index)
    value = format["component"]("dof_values", format_s)
    lines_s = [format_assign(name, value)]
    lines_r += generate_loop(lines_s, loop_vars_s, Indent, format)
    code += generate_loop(lines_r, loop_vars_r, Indent, format)

    code += ["", format["comment"]("Delete pointer.")]
    code += [Indent.indent(format["delete pointer"]("dof_values", ""))]

    # Generate bode (no need to remove unused)
    return "\n".join(code)

def _evaluate_basis_derivatives(data_list):
    """Evaluate the derivatives of an element basisfunction at a point. The values are
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
    Nedelec (broken?)"""

    #set_format(format)

    # Init return code and indent object
    code = []
    Indent = IndentControl()

    # Get the element cell domain, geometric and topological dimension.
    element_cell_domain = data_list[0]["cell_domain"]
    geometric_dimension = data_list[0]["geometric_dimension"]
    topological_dimension = data_list[0]["topological_dimension"]

    # Get code snippets for Jacobian, Inverse of Jacobian and mapping of
    # coordinates from physical element to the FIAT reference element.
    # FIXME: KBO: Change this when supporting R^2 in R^3 elements.
    code += [Indent.indent(format["jacobian and inverse"](geometric_dimension))]
    code += ["", Indent.indent(format["fiat coordinate map"](element_cell_domain))]

    # Compute number of derivatives that has to be computed, and declare an array to hold
    # the values of the derivatives on the reference element
    code += [""]
    code += _compute_num_derivatives(topological_dimension, Indent, format)

    # Generate all possible combinations of derivatives
    code += _generate_combinations(topological_dimension, Indent, format)

    # Generate the transformation matrix
    code += _generate_transform(element_cell_domain, Indent, format)

    # Reset all values
    code += _reset_values(data_list, Indent, format)

    if len(data_list) == 1:
        data = data_list[0]

        # Map degree of freedom to local degree.
        code += [_map_dof(0, Indent, format)]

        # Generate element code.
        code += _generate_element_code(data, 0, Indent, format)

    # If the element is of type MixedElement (including Vector- and TensorElement).
    else:
        # Generate element code, for all sub-elements.
        code += _mixed_elements(data_list, Indent, format)
    code = remove_unused("\n".join(code))
    return code

def _compute_num_derivatives(topological_dimension, Indent, format):
    "Computes the number of derivatives of order 'n' as: element.cell_shape()^n."

    code = [format["comment"]("Compute number of derivatives.")]
    code += [format["declaration"](Indent.indent(format["uint declaration"]), format["num derivatives"],\
              format["floating point"](1))]

    loop_vars = [(format["free indices"][0], 0, format["argument derivative order"])]
    lines = [format["imul"](format["num derivatives"], format["floating point"](topological_dimension))]

    code += generate_loop(lines, loop_vars, Indent, format)

    return code

def _generate_combinations(topological_dimension, Indent, format):
    "Generate all possible combinations of derivatives of order 'n'"

    # Use code from format
    code = ["", Indent.indent(format["combinations"]\
            % {"combinations": format["derivative combinations"],\
               "topological_dimension-1": topological_dimension-1,\
               "num_derivatives" : format["num derivatives"],\
               "n": format["argument derivative order"]})]
    return code

def _generate_transform(element_cell_domain, Indent, format):
    """Generate the transformation matrix, whic is used to transform derivatives from reference
    element back to the physical element."""

    # Generate code to construct the inverse of the Jacobian
    if (element_cell_domain in ["interval", "triangle", "tetrahedron"]):
        code = ["", Indent.indent(format["transform snippet"][element_cell_domain]\
        % {"transform": format["transform matrix"],\
           "num_derivatives" : format["num derivatives"],\
           "n": format["argument derivative order"],\
           "combinations": format["derivative combinations"],\
           "K":format["transform Jinv"]})]
    else:
        error("Cannot generate transform for shape: %s" % element_cell_domain)

    return code

def _reset_values(data_list, Indent, format):
    "Reset all components of the 'values' array as it is a pointer to an array."
    format_assign = format["assign"]
    code = ["", Indent.indent(format["comment"]("Reset values. Assuming that values is always an array."))]

    # Get value shape and reset values. This should also work for TensorElement,
    # scalar are empty tuples, therefore (1,) in which case value_shape = 1.
    value_shape = sum(sum(data["value_shape"] or (1,)) for data in data_list)

    # Only multiply by value shape if different from 1.
    if value_shape == 1:
        num_vals = format["num derivatives"]
    else:
        num_vals = format["multiply"]([format["floating point"](value_shape), format["num derivatives"]])
    name = format["component"](format["argument values"], format["free indices"][0])
    loop_vars = [(format["free indices"][0], 0, num_vals)]
    lines = [format_assign(name, format["floating point"](0))]
    code += generate_loop(lines, loop_vars, Indent, format)

    return code + [""]

def _generate_element_code(data, sum_value_dim, Indent, format):
    "Generate code for each basis element"

    code = []

    # Compute basisvalues, from evaluatebasis.py
    code += _compute_basisvalues(data, Indent, format)

    # Tabulate coefficients
    code += _tabulate_coefficients(data, Indent, format)

    # Tabulate coefficients for derivatives
    code += _tabulate_dmats(data, Indent, format)

    # Compute the derivatives of the basisfunctions on the reference (FIAT) element,
    # as the dot product of the new coefficients and basisvalues
    code += _compute_reference_derivatives(data, Indent, format)

    # Transform derivatives to physical element by multiplication with the transformation matrix
    code += _transform_derivatives(data, sum_value_dim, Indent, format)

    # Delete pointers
    code += _delete_pointers(data, Indent, format)

    return code

def _mixed_elements(data_list, Indent, format):
    "Generate code for each sub-element in the event of vector valued elements or mixed elements"

    # Prefetch formats to speed up code generation
    format_dof_map_if = format["dof map if"]
    format_if         = format["if"]

    sum_value_dim = 0
    sum_space_dim = 0

    # Init return code.
    code = []

    # Generate code for each element
    for data in data_list:

        # Get value and space dimension (should be tensor ready).
        value_dim = sum(data["value_shape"] or (1,))
        space_dim = data["space_dimension"]

        # Generate map from global to local dof
        element_code = [_map_dof(sum_space_dim, Indent, format)]

        # Generate code for basis element
        element_code += _generate_element_code(data, sum_value_dim, Indent, format)

        # Increase indentation, indent code and decrease indentation.
        Indent.increase()
        if_code = remove_unused(Indent.indent("\n".join(element_code)))
        # if_code = Indent.indent("\n".join(element_code))
        Indent.decrease()

        # Create if statement and add to code.
        code += [format_if(format_dof_map_if(sum_space_dim, sum_space_dim + space_dim -1), if_code)]

        # Increase sum of value dimension, and space dimension
        sum_value_dim += value_dim
        sum_space_dim += space_dim

    return code

def _tabulate_dmats(data, Indent, format):
    "Tabulate the derivatives of the polynomial base"

    code = []

    # Prefetch formats to speed up code generation
    format_table          = format["static const float declaration"]
    format_dmats          = format["dmats"]
    format_component  = format["component"]
    format_assign = format["assign"]

    # Get derivative matrices (coefficients) of basis functions, computed by FIAT at compile time
    derivative_matrices = data["dmats"]

    code += [Indent.indent(format["comment"]("Tables of derivatives of the polynomial base (transpose)."))]

    # Generate tables for each spatial direction
    for i, dmat in enumerate(derivative_matrices):

        # Extract derivatives for current direction (take transpose, FIAT_NEW PolynomialSet.tabulate())
        matrix = numpy.transpose(dmat)

        # Get shape and check dimension (This is probably not needed)
        shape = numpy.shape(matrix)
        ffc_assert(shape[0] == shape[1] == data["num_expansion_members"], "Something is wrong with the shape of dmats.")

        # Declare varable name for coefficients
        name = format_component(format_table + format_dmats(i), [shape[0], shape[1]])
        value = tabulate_matrix(matrix, format)
        code += [format_assign(Indent.indent(name), Indent.indent(value)), ""]

    return code

def _reset_dmats(shape_dmats, indices, Indent, format):
    format_assign = format["assign"]
    code = [format["comment"]("Resetting dmats values to compute next derivative.")]

    loop_vars = [(indices[0], 0, shape_dmats[0]), (indices[1], 0, shape_dmats[1])]
    dmats_old = format["component"](format["dmats"](""), [indices[0], indices[1]])
    lines = [format_assign(dmats_old, format["floating point"](0.0))]
    lines += [format["if"](indices[0] + format["is equal"] + indices[1],\
              format["assign"](Indent.indent(dmats_old), format["floating point"](1.0)))]
    code += generate_loop(lines, loop_vars, Indent, format)
    return code

def _update_dmats(shape_dmats, indices, Indent, format):
    format_assign = format["assign"]
    code = [format["comment"]("Updating dmats_old with new values and resetting dmats.")]
    dmats = format["component"](format["dmats"](""), [indices[0], indices[1]])
    dmats_old = format["component"](format["dmats old"], [indices[0], indices[1]])
    loop_vars = [(indices[0], 0, shape_dmats[0]), (indices[1], 0, shape_dmats[1])]
    lines = [format_assign(dmats_old, dmats), format_assign(dmats, format["floating point"](0.0))]
    code += generate_loop(lines, loop_vars, Indent, format)
    return code

def _compute_dmats(num_dmats, shape_dmats, available_indices, Indent, format):

    s, t, u = available_indices

    # Reset dmats_old
    code = _reset_dmats(shape_dmats, [t, u], Indent, format)
    code += ["", format["comment"]("Looping derivative order to generate dmats.")]

    # Set dmats matrix equal to dmats_old
    lines = _update_dmats(shape_dmats, [t, u], Indent, format)
    loop_vars = [(s, 0, format["argument derivative order"])]

    lines += ["", format["comment"]("Update dmats using an inner product.")]
    # Create dmats matrix by multiplication
    for i in range(num_dmats):
        lines += _dmats_product(shape_dmats, s, i, [t, u], Indent, format)

    code += generate_loop(lines, loop_vars, Indent, format)

    return code

def _dmats_product(shape_dmats, index, i, indices, Indent, format):

    t, u = indices
    tu = t + u
    loop_vars = [(t, 0, shape_dmats[0]), (u, 0, shape_dmats[1])]
    dmats = format["component"](format["dmats"](""), [t, u])
    dmats_old = format["component"](format["dmats old"], [tu, u])
    value = format["multiply"]([format["component"](format["dmats"](i), [t, tu]), dmats_old])
    name = Indent.indent(format["iadd"](dmats, value))
    lines = generate_loop([name], [(tu, 0, shape_dmats[0])], Indent, format)

    code = [format["if"](index + format["is equal"] + str(i),\
            "\n".join(generate_loop(lines, loop_vars, Indent, format)))]

    return code

def _compute_component(shape_dmats, index, component, indices, Indent, format):
    loop_vars = [(indices[0], 0, shape_dmats[0]),(indices[1], 0, shape_dmats[1])]

    if component == 0:
        access = index
    elif component == 1:
        access = format["add"]([format["num derivatives"], index])
    else:
        mul = format["multiply"]([format["floating point"](component), format["num derivatives"]])
        access = format["add"]([mul, index])

    name = format["component"](format["reference derivatives"], access)
    coeffs = format["component"](format["coefficients"](component), [format["local dof"], indices[0]])
    dmats = format["component"](format["dmats"](""), [indices[0], indices[1]])
    basis = format["component"](format["basisvalues"], indices[1])
    value = format["multiply"]([coeffs, dmats, basis])
    return generate_loop([format["iadd"](name, value)], loop_vars, Indent, format)

def _compute_reference_derivatives(data, Indent, format):
    """Compute derivatives on the reference element by recursively multiply coefficients with
    the relevant derivatives of the polynomial base until the requested order of derivatives
    has been reached. After this take the dot product with the basisvalues."""

    code = []

    # Prefetch formats to speed up code generation
    format_comment          = format["comment"]
    format_float            = format["float declaration"]
#    format_coeff            = format["coefficient scalar"]
#    format_new_coeff        = format["new coefficient scalar"]
#    format_secondary_index  = format["secondary index"]
#    format_floating_point   = format["floating point"]
#    format_num_derivatives  = format["num derivatives"]
#    format_loop             = format["loop"]
#    format_block_begin      = format["block begin"]
#    format_block_end        = format["block end"]
#    format_coefficients     = format["coefficients"]
#    format_dof              = format["local dof"]
#    format_n                = format["argument derivative order"]
#    format_derivatives      = format["reference derivatives"]
    format_component    = format["component"]
#    format_add              = format["add"]
#    format_multiply         = format["multiply"]
#    format_inv              = format["inverse"]
#    format_det              = format["determinant"]
#    format_group            = format["grouping"]
#    format_basisvalue       = format["basisvalue"]
    format_tmp              = format["tmp ref value"]
#    format_tmp_access       = format["tmp access"]
#    format_table            = format["table declaration"]
    format_dmats            = format["dmats"]
    format_assign = format["assign"]

    format_r, format_s, format_t, format_u = format["free indices"]

    # Get number of components, change for tensor valued elements.
    shape = data["value_shape"]
    if shape == ():
        num_components = 1
    elif len(shape) == 1:
        num_components = shape[0]
    else:
        error("Tensor valued elements are not supported yet: %s" % data["family"])

    # Get shape of derivative matrix (they should all have the same shape).
    shape_dmats = numpy.shape(data["dmats"][0])

    code += [Indent.indent(format_comment("Compute reference derivatives"))]

    # Declare pointer to array that holds derivatives on the FIAT element
    code += [Indent.indent(format_comment("Declare pointer to array of derivatives on FIAT element"))]
    # The size of the array of reference derivatives is equal to the number of derivatives
    # times the number of components of the basis element
    if (num_components == 1):
        num_vals = format["num derivatives"]
    else:
        num_vals = format["multiply"]([format["floating point"](num_components), format["num derivatives"]])
    code += [format_assign(Indent.indent(format_float + format["pointer"] + format["reference derivatives"]),\
              format["component"](format["new"] + format_float, num_vals))]
    # Reset values of reference derivatives.
    name = format["component"](format["reference derivatives"], format_r)
    lines = [format_assign(name, format["floating point"](0))]
    code += generate_loop(lines, [(format_r, 0, num_vals)], Indent, format)

    code += [""]

    # Declare matrix of dmats (which will hold the matrix product of all combinations)
    # and dmats_old which is needed in order to perform the matrix product.
    code += [Indent.indent(format_comment("Declare derivative matrix (of polynomial basis)."))]
    matrix = numpy.eye(shape_dmats[0])
    name = format_component(format_float + format_dmats(""), [shape_dmats[0], shape_dmats[1]])
    value = tabulate_matrix(matrix, format)
    code += [format_assign(Indent.indent(name), Indent.indent(value)), ""]
    code += [Indent.indent(format_comment("Declare (auxiliary) derivative matrix (of polynomial basis)."))]
    name = format_component(format_float + format["dmats old"], [shape_dmats[0], shape_dmats[1]])
    code += [format_assign(Indent.indent(name), Indent.indent(value)), ""]

    # Loop all derivatives and compute value of the derivative as:
    # deriv_on_ref[r] = coeff[dof][s]*dmat[s][t]*basis[t]
    code += [Indent.indent(format_comment("Loop possible derivatives."))]
    loop_vars = [(format_r, 0, format["num derivatives"])]
    # Compute dmats as a recursive matrix product
    lines = _compute_dmats(len(data["dmats"]), shape_dmats, [format_s, format_t, format_u], Indent, format)
    # Compute derivatives for all components
    for i in range(num_components):
        lines += _compute_component(shape_dmats, format_r, i, [format_s, format_t], Indent, format)

    code += generate_loop(lines, loop_vars, Indent, format)

    # Apply transformation if applicable.
    mapping = data["mapping"]
    if mapping == "affine":
        pass
    elif mapping == "contravariant piola":
        code += ["", Indent.indent(format["comment"]\
                ("Using contravariant Piola transform to map values back to the physical element"))]
        # Get temporary values before mapping.
        code += [format["const float declaration"](Indent.indent(format_tmp(i)),\
                  format_component(format["reference derivatives"], i)) for i in range(num_components)]

        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [format_tmp(j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create Jacobian.
            jacobian_row = [format["transform"]("J", j, i, None) for j in range(topological_dimension)]

            # Create inner product and multiply by inverse of Jacobian.
            inner = [format["multiply"]([jacobian_row[j], basis_col[j]]) for j in range(topological_dimension)]
            sum_ = format["grouping"](format["add"](inner))
            value = format["multiply"]([format["inverse"](format["det(J)"]("")), sum_])
            name = format_component(format["reference derivatives"], i)
            code += [format_assign(name, value)]
    elif mapping == "covariant piola":
        code += ["", Indent.indent(format["comment"]\
                ("Using covariant Piola transform to map values back to the physical element"))]
        # Get temporary values before mapping.
        code += [format["const float declaration"](Indent.indent(format_tmp(i)),\
                  format_component(format["reference derivatives"], i)) for i in range(num_components)]
        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [format_tmp(j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create inverse of Jacobian.
            inv_jacobian_row = [format["transform"]("JINV", j, i, None) for j in range(topological_dimension)]

            # Create inner product of basis values and inverse of Jacobian.
            inner = [format["multiply"]([inv_jacobian_row[j], basis_col[j]]) for j in range(topological_dimension)]
            value = format["grouping"](format["add"](inner))
            name = format_component(format["reference derivatives"], i)
            code += [format_assign(name, value)]
    else:
        error("Unknown mapping: %s" % mapping)

    return code + [""]

def _transform_derivatives(data, sum_value_dim, Indent, format):
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
    format_component    = format["component"]
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

    # Get number of components, change for tensor valued elements.
    shape = data["value_shape"]
    if shape == ():
        num_components = 1
    elif len(shape) == 1:
        num_components = shape[0]
    else:
        error("Tensor valued elements are not supported yet: %s" % data["family"])

    # Compute offset in array values if any
    for i in range(num_components):
        if (sum_value_dim + i == 0):
            name = format_component(format_values, "row")
        elif (sum_value_dim + i == 1):
            name = format_component(format_values, format_add([format_num_derivatives, "row"]))
        else:
            offset_name = format_multiply(["%d" %(sum_value_dim + i), format_num_derivatives])
            name = format_component(format_values, format_add([offset_name, "row"]))
        if (i == 0):
            value = format_multiply([format_component(format_transform, ["row","col"]),\
                    format_component(format_derivatives, "col")])
        elif (i == 1):
            value = format_multiply([format_component(format_transform, ["row","col"]),\
                    format_component(format_derivatives, format_add([format_num_derivatives, "col"]))])
        else:
            offset_value = format_multiply(["%d" %i, format_num_derivatives])

            value = format_multiply([format_component(format_transform, ["row","col"]),\
                    format_component(format_derivatives, format_add([offset_value, "col"]))])

        # Compute values
        code += [Indent.indent(format["iadd"](name, value))]

    # Decrease indentation
    Indent.decrease()
    code += [Indent.indent(format["block end"])]

    # Decrease indentation
    Indent.decrease()

    code += [Indent.indent(format["block end"])]

    return code

def _delete_pointers(data, Indent, format):
    "Delete the pointers to arrays."

    code = []
    format_r = format["free indices"][0]


    # Delete pointers
    code += ["", Indent.indent(format["comment"]("Delete pointer to array of derivatives on FIAT element"))]
    code += [Indent.indent(format["delete pointer"](format["reference derivatives"], "")), ""]

    code += [Indent.indent(format["comment"]("Delete pointer to array of combinations of derivatives and transform"))]
    loop_vars = [(format_r, 0, format["num derivatives"])]
    lines =  [format["delete pointer"](format["derivative combinations"], format["component"]("", format_r))]
    lines += [format["delete pointer"](format["transform matrix"], format["component"]("", format_r))]
    code += generate_loop(lines, loop_vars, Indent, format)

    code += [Indent.indent(format["delete pointer"](format["derivative combinations"], ""))]
    code += [Indent.indent(format["delete pointer"](format["transform matrix"], ""))]

    return code

