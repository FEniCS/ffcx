"""Code generation for evaluation of derivatives of finite element basis values.
This module generates code which is more or less a C++ representation of the code
found in FIAT_NEW."""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-04-16"
__copyright__ = "Copyright (C) 2007-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-08

# Python modules
import math
import numpy

# FFC modules
from ffc.log import error, ffc_assert
from ffc.evaluatebasis import _map_dof
from ffc.evaluatebasis import _compute_basisvalues
from ffc.evaluatebasis import _tabulate_coefficients
from ffc.cpp import tabulate_matrix, IndentControl, remove_unused
from ffc.quadrature.quadraturegenerator_utils import generate_loop
#from ffc.quadrature.symbolics import set_format

# Temporary import
from cpp import format_old as format

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

    # Get the element cell domain and check if it is the same for all elements.
    element_cell_domain = data_list[0]["cell_domain"]
    # FIXME: KBO: This should already be checked elsewhere.
    ffc_assert(all(element_cell_domain == data["cell_domain"] for data in data_list),\
               "The element cell domain must be the same for all sub elements: " + repr(data_list))

    # Create mapping from physical element to the FIAT reference element
    # (from codesnippets.py).
    # FIXME: KBO: Change this when supporting R^2 in R^3 elements.
    code += [Indent.indent(format["coordinate map FIAT"](element_cell_domain))]

    # Get the topological dimension.
    # FIXME: If the elements are defined on the same cell domain they also have
    # the same topological dimension so there is no need to check this.
    topological_dimension = data_list[0]["topological_dimension"]
    ffc_assert(all(topological_dimension == data["topological_dimension"] for data in data_list),\
               "The topological dimension must be the same for all sub elements: " + repr(data_list))

    # Compute number of derivatives that has to be computed, and declare an array to hold
    # the values of the derivatives on the reference element
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
        code += _map_dof(0, Indent, format)

        # Generate element code.
        code += _generate_element_code(data, 0, Indent, format)

    # If the element is of type MixedElement (including Vector- and TensorElement).
    else:
        # Generate element code, for all sub-elements.
        code += _mixed_elements(data_list, Indent, format)

    lines = format["generate body"](code)
#    code = remove_unused(lines)
    return lines
    return code

def _compute_num_derivatives(topological_dimension, Indent, format):
    "Computes the number of derivatives of order 'n' as: element.cell_shape()^n."

    code = ["", format["comment"]("Compute number of derivatives.")]
    code += [(Indent.indent(format["uint declaration"] + format["num derivatives"]),\
              format["floating point"](1))]

    loop_vars = [(format["free secondary indices"][0], 0, format["argument derivative order"])]
    lines = [format["times equal"](format["num derivatives"],\
             format["floating point"](topological_dimension))]

    code += generate_loop(lines, loop_vars, Indent, format)

    return code

def _generate_combinations(topological_dimension, Indent, format):
    "Generate all possible combinations of derivatives of order 'n'"

    # Use code from codesnippets.py
    code = ["", Indent.indent(format["snippet combinations"])\
            % {"combinations": format["derivative combinations"],\
               "topological_dimension-1": topological_dimension-1,\
               "num_derivatives" : format["num derivatives"],\
               "n": format["argument derivative order"]}]
    return code

def _generate_transform(element_cell_domain, Indent, format):
    """Generate the transformation matrix, whic is used to transform derivatives from reference
    element back to the physical element."""

    # Generate code to construct the inverse of the Jacobian, use code from codesnippets.py
    if (element_cell_domain in ["interval", "triangle", "tetrahedron"]):
        code = ["", Indent.indent(format["snippet transform"](element_cell_domain))\
        % {"transform": format["transform matrix"],\
           "num_derivatives" : format["num derivatives"],\
           "n": format["argument derivative order"],\
           "combinations": format["derivative combinations"],\
           "K":format["transform Jinv"]}]
    else:
        error("Cannot generate transform for shape: %s" % element_cell_domain)

    return code

def _reset_values(data_list, Indent, format):
    "Reset all components of the 'values' array as it is a pointer to an array."

    code = ["", Indent.indent(format["comment"]("Reset values. Assuming that values is always an array."))]

    # Get value shape and reset values. This should also work for TensorElement,
    # scalar are empty tuples, therefore (1,) in which case value_shape = 1.
    value_shape = sum(sum(data["value_shape"] or (1,)) for data in data_list)

    # Only multiply by value shape if different from 1.
    if value_shape == 1:
        num_vals = format["num derivatives"]
    else:
        num_vals = format["multiply"]([format["floating point"](value_shape), format["num derivatives"]])
    name = format["argument values"] + format["array access"](format["free secondary indices"][0])
    loop_vars = [(format["free secondary indices"][0], 0, num_vals)]
    lines = [(name, format["floating point"](0))]
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

def _mixed_elements(data, Indent, format):
    "Generate code for each sub-element in the event of vector valued elements or mixed elements"

    code = []

    # Prefetch formats to speed up code generation
    format_block_begin  = format["block begin"]
    format_block_end    = format["block end"]
    format_dof_map_if   = format["dof map if"]

    # Extract basis elements, and determine number of elements
    elements = data.extract_elements()
    num_elements = len(datas)

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
        code += _dof_map(sum_space_dim, Indent, format)

        # Generate code for basis element
        code += _generate_element_code(basis_element, sum_value_dim, Indent, format)

        # Decrease indentation, finish block - end element code
        Indent.decrease()
        code += [Indent.indent(format_block_end)] + [""]

        # Increase sum of value dimension, and space dimension
        sum_value_dim += value_dim
        sum_space_dim += space_dim

    return code

def _tabulate_dmats(data, Indent, format):
    "Tabulate the derivatives of the polynomial base"

    code = []

    # Prefetch formats to speed up code generation
    format_table          = format["table declaration"]
    format_dmats          = format["dmats table"]
    format_matrix_access  = format["matrix access"]

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
        name = format_table + format_dmats(i) + format_matrix_access(shape[0], shape[1])
        value = tabulate_matrix(matrix, format)
        code += [(Indent.indent(name), Indent.indent(value)), ""]

    return code

def _reset_dmats(shape_dmats, indices, Indent, format):
    code = [format["comment"]("Resetting dmats values to compute next derivative.")]

    loop_vars = [(indices[0], 0, shape_dmats[0]), (indices[1], 0, shape_dmats[1])]
    dmats_old = format["dmats"] + format["matrix access"](indices[0], indices[1])
    lines = [(dmats_old, format["floating point"](0.0))]
    lines += [format["if"] + format["grouping"](indices[0] + format["is equal"] + indices[1])]
    Indent.increase()
    lines += [(Indent.indent(dmats_old), format["floating point"](1.0))]
    Indent.decrease()
    code += generate_loop(lines, loop_vars, Indent, format)
    return code

def _update_dmats(shape_dmats, indices, Indent, format):
    code = [format["comment"]("Updating dmats_old with new values and resetting dmats.")]
    dmats = format["dmats"] + format["matrix access"](indices[0], indices[1])
    dmats_old = format["dmats old"] + format["matrix access"](indices[0], indices[1])
    loop_vars = [(indices[0], 0, shape_dmats[0]), (indices[1], 0, shape_dmats[1])]
    lines = [(dmats_old, dmats), (dmats, format["floating point"](0.0))]
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

    code = [format["if"] + format["grouping"](index + format["is equal"] + str(i))]
    code += [format["block begin"]]
    Indent.increase()
    loop_vars = [(t, 0, shape_dmats[0]), (u, 0, shape_dmats[1])]
    dmats = format["dmats"] + format["matrix access"](t, u)
    dmats_old = format["dmats old"] + format["matrix access"](tu, u)
    value = format["multiply"]([format["dmats table"](i) + format["matrix access"](t, tu), dmats_old])
    name = Indent.indent(format["add equal"](dmats, value))
    lines = generate_loop([name], [(tu, 0, shape_dmats[0])], Indent, format)
    code += generate_loop(lines, loop_vars, Indent, format)
    Indent.decrease()
    code += [format["block end"]]

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

    name = format["reference derivatives"] + format["array access"](access)
    coeffs = format["coefficients table"](component) + format["matrix access"](format["local dof"], indices[0])
    dmats = format["dmats"] + format["matrix access"](indices[0], indices[1])
    basis = format["basisvalues"](indices[1])
    value = format["multiply"]([coeffs, dmats, basis])
    return generate_loop([format["add equal"](name, value)], loop_vars, Indent, format)

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
#    format_coefficients     = format["coefficients table"]
#    format_dof              = format["local dof"]
#    format_n                = format["argument derivative order"]
#    format_derivatives      = format["reference derivatives"]
    format_matrix_access    = format["matrix access"]
    format_array_access     = format["array access"]
#    format_add              = format["add"]
#    format_multiply         = format["multiply"]
#    format_inv              = format["inverse"]
#    format_det              = format["determinant"]
#    format_group            = format["grouping"]
#    format_basisvalue       = format["basisvalue"]
    format_tmp              = format["tmp declaration"]
    format_tmp_access       = format["tmp access"]
#    format_table            = format["table declaration"]
    format_dmats            = format["dmats"]

    format_r, format_s, format_t, format_u = format["free secondary indices"]

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
    code += [(Indent.indent(format_float + format["pointer"] + format["reference derivatives"]),\
              format["new"] + format_float + format["array access"](num_vals))]
    # Reset values of reference derivatives.
    name = format["reference derivatives"] + format["array access"](format_r)
    lines = [(name, format["floating point"](0))]
    code += generate_loop(lines, [(format_r, 0, num_vals)], Indent, format)

    code += [""]

    # Declare matrix of dmats (which will hold the matrix product of all combinations)
    # and dmats_old which is needed in order to perform the matrix product.
    code += [Indent.indent(format_comment("Declare derivative matrix (of polynomial basis)."))]
    matrix = numpy.eye(shape_dmats[0])
    name = format_float + format_dmats + format_matrix_access(shape_dmats[0], shape_dmats[1])
    value = tabulate_matrix(matrix, format)
    code += [(Indent.indent(name), Indent.indent(value)), ""]
    code += [Indent.indent(format_comment("Declare (auxiliary) derivative matrix (of polynomial basis)."))]
    name = format_float + format["dmats old"] + format_matrix_access(shape_dmats[0], shape_dmats[1])
    code += [(Indent.indent(name), Indent.indent(value)), ""]

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
        code += [(Indent.indent(format_tmp(0, i)),\
                  format["reference derivatives"] + format_array_access(i)) for i in range(num_components)]

        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [format_tmp_access(0, j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create Jacobian.
            jacobian_row = [format["transform"]("J", j, i, None) for j in range(topological_dimension)]

            # Create inner product and multiply by inverse of Jacobian.
            inner = [format["multiply"]([jacobian_row[j], basis_col[j]]) for j in range(topological_dimension)]
            sum_ = format["grouping"](format["add"](inner))
            value = format["multiply"]([format["inverse"](format["determinant"](None)), sum_])
            name = format["reference derivatives"] + format_array_access(i)
            code += [(name, value)]
    elif mapping == "covariant piola":
        code += ["", Indent.indent(format["comment"]\
                ("Using covariant Piola transform to map values back to the physical element"))]
        # Get temporary values before mapping.
        code += [(Indent.indent(format_tmp(0, i)),\
                  format["reference derivatives"] + format_array_access(i)) for i in range(num_components)]
        # Create names for inner product.
        topological_dimension = data["topological_dimension"]
        basis_col = [format_tmp_access(0, j) for j in range(topological_dimension)]
        for i in range(num_components):
            # Create inverse of Jacobian.
            inv_jacobian_row = [format["transform"]("JINV", j, i, None) for j in range(topological_dimension)]

            # Create inner product of basis values and inverse of Jacobian.
            inner = [format["multiply"]([inv_jacobian_row[j], basis_col[j]]) for j in range(topological_dimension)]
            value = format["grouping"](format["add"](inner))
            name = format["reference derivatives"] + format_array_access(i)
            code += [(name, value)]
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

def _delete_pointers(data, Indent, format):
    "Delete the pointers to arrays."

    code = []

    # Delete pointers
    code += [Indent.indent(format["comment"]("Delete pointer to array of derivatives on FIAT element"))]
    code += [Indent.indent(format["delete"] + format["array access"]("") + " " +\
                           format["reference derivatives"] + format["end line"])] + [""]

    code += [Indent.indent(format["comment"]("Delete pointer to array of combinations of derivatives and transform"))]


    code += [Indent.indent(format["loop"]("row", 0, format["num derivatives"]))]
    code += [Indent.indent(format["block begin"])]
    # Increase indentation
    Indent.increase()

    code += [Indent.indent(format["delete"] + format["array access"]("") + " " +\
                           format["derivative combinations"] + format["array access"]("row") + format["end line"])]

    code += [Indent.indent(format["delete"] + format["array access"]("") + " " +\
                           format["transform matrix"] + format["array access"]("row") + format["end line"])]

    # Decrease indentation
    Indent.decrease()
    code += [Indent.indent(format["block end"])] + [""]

    code += [Indent.indent(format["delete"] + format["array access"]("") + " " +\
                           format["derivative combinations"] + format["end line"])]

    code += [Indent.indent(format["delete"] + format["array access"]("") + " " +\
                           format["transform matrix"] + format["end line"])]

    return code

