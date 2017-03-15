# -*- coding: utf-8 -*-
"""Work in progress translation of FFC evaluatebasis code to uflacs CNodes format."""

from six import string_types
import math
import numpy

from ffc.log import error
from ffc.uflacs.backends.ufc.utils import generate_error

from ffc.uflacs.backends.ufc.evaluatebasis import generate_expansion_coefficients, generate_compute_basisvalues


# Used for various arrays in this file
int_type = "int"


def generate_evaluate_reference_basis_derivatives(L, data, parameters):
    # Cutoff for feature to disable generation of this code (consider removing after benchmarking final result)
    if isinstance(data, string_types):
        msg = "evaluate_reference_basis_derivatives: %s" % (data,)
        return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

    # Get some known dimensions
    element_cellname = data["cellname"]
    tdim = data["topological_dimension"]
    max_degree = data["max_degree"]
    reference_value_size = data["reference_value_size"]
    num_dofs = len(data["dofs_data"])

    # Output argument
    reference_values = L.Symbol("reference_values")

    # Input arguments
    order = L.Symbol("order")
    num_points = L.Symbol("num_points")
    X = L.Symbol("X")

    ref_values = L.FlattenedArray(reference_values, dims=(num_points, num_dofs, reference_value_size))  # FIXME: Dimensions?

    # Loop indices
    ip = L.Symbol("ip")
    k = L.Symbol("k")
    c = L.Symbol("c")
    r = L.Symbol("r")

    # Cutoff to evaluate_basis for order 0
    call_eval_basis_code = [
        L.Comment("Call evaluate_basis if order of derivatives is equal to zero"),
        L.If(L.EQ(order, 0), [
            L.Call("evaluate_reference_basis", (reference_values, num_points, X)),
            L.Return()
            ])
        ]

    # Tabulate dmats tables for all dofs and all derivative directions
    dmats_code, dmats_names = generate_tabulate_dmats(L, data["dofs_data"])

    # Generate code with static tables of expansion coefficients
    tables_code, coefficients_for_dof = generate_expansion_coefficients(L, data["dofs_data"])

    # Precompute number of derivatives for each order
    all_num_derivatives = L.Symbol("all_num_derivatives")
    count_derivs_code = [
        L.Comment("Number of derivatives for each order"),
        L.ArrayDecl("static const " + int_type, all_num_derivatives, max_degree+1,
                    values=[tdim**n for n in range(max_degree+1)]),
        ]
    num_derivatives = all_num_derivatives[order]
    # FIXME: num_derivatives was computed runtime for a reason: to support zero reset for order > max_degree

    reset_values_code = [
        L.Comment("Reset reference_values[:] to 0"),
        L.ForRange(k, 0, num_points * num_dofs * num_derivatives * reference_value_size, body=
            L.Assign(reference_values[k], 0.0))
        ]

    order_cutoff_code = [
        L.Comment("If order of derivatives is greater than the maximum polynomial degree, return zeros"),
        L.If(L.GT(order, max_degree),
             L.Return()),
        ]

    # Collect pieces
    setup_code = (
        dmats_code
        + tables_code
        + call_eval_basis_code
        + count_derivs_code
        + reset_values_code
        + order_cutoff_code
        )

    # If max_degree is zero, we're done here
    if max_degree == 0:
        return setup_code

    # Generate code to compute tables of basisvalues
    basisvalues_code, basisvalues_for_degree, need_fiat_coordinates = \
        generate_compute_basisvalues(L, data["dofs_data"], element_cellname, tdim, X, ip)


    # TODO: Access dmats as dmats_names[idof][ideriv]

    # Accumulate products of basisvalues and coefficients into values
    accumulation_code = generate_accumulation_code(L,
        data["dofs_data"],
        num_derivatives,
        coefficients_for_dof)


    # Generate geo code.
    geometry_code = [] # _geometry_related_code(L, data, tdim, gdim, element_cellname)  # FIXME

    # Generate all possible combinations of derivatives.
    # FIXME: What is this really doing? Translation below.
    #        Probably need to look at code generation for applying mappings that was deleted.
    combinations_code = _generate_combinations(L, tdim, "", max_degree)

    # Create code for all basis values (dofs).
    dof_cases = []
    for dof_data in data["dofs_data"]:
        case = _generate_dof_code(data, dof_data)
        dof_cases.append(case)

    # Define loop over points
    final_loop_code = [
            L.ForRange(ip, 0, num_points,
                   body=basisvalues_code + accumulation_code)
        ]

    # Stitch it all together
    code = (
        setup_code
        + geometry_code,
        + combinations_code,
        + final_loop_code
        )
    return code


def _generate_dof_code(data, dof_data):
    "Generate code for a basis."
    code = []

    # Compute the derivatives of the basisfunctions on the reference (FIAT) element,
    # as the dot product of the new coefficients and basisvalues.
    code += _compute_reference_derivatives(data, dof_data) # FIXME: Convert this

    # Transform derivatives to physical element by multiplication with the transformation matrix.
    code += _transform_derivatives(data, dof_data)  # TODO: Extract to separate ufc function

    code += [format["switch"](format["argument basis num"], dof_cases)]  # TODO: Loop instead of switch

    return code


def generate_tabulate_dmats(L, dofs_data):
    "Tabulate the derivatives of the polynomial base"

    alignas = 32

    # TODO: Do we want to share some dmats between dofs?
    # Or make dofs a dimension in the tables instead of multiple tables?
    # Need to look at data to determine, currently just transforming code blindly.

    code = [L.Comment("Tables of derivatives of the polynomial base (transpose).")]

    dmats_names = []
    #import ipdb; ipdb.set_trace()
    for idof, dof_data in enumerate(dofs_data):
        # Get derivative matrices (coefficients) of basis functions, computed by FIAT at compile time.
        derivative_matrices = dof_data["dmats"]

        dmats_names.append([])

        # Generate tables for each spatial direction.
        for i, dmat in enumerate(derivative_matrices):

            # Extract derivatives for current direction
            # (take transpose, FIAT_NEW PolynomialSet.tabulate()).
            matrix = numpy.transpose(dmat)

            # Get shape and check dimension (This is probably not needed).
            shape = numpy.shape(matrix)
            if not (shape[0] == shape[1] == dof_data["num_expansion_members"]):
                error("Something is wrong with the shape of dmats.")

            # Declare varable name for coefficients.
            name = L.Symbol("dmats%d_%d" % (i, idof))
            code += [L.ArrayDecl("static const double", name, shape, values=matrix, alignas=alignas)]
            dmats_names[-1].append(name)

    return code, dmats_names


def _generate_combinations(L, tdim, dummy, max_degree):
    # don't remember what the "" argument, here dummy, was supposed to be

    # Input argument
    order = L.Symbol("order")

    num_derivatives = L.Symbol("num_derivatives") # FIXME: To be computed from order

    max_num_derivatives = tdim**max_degree

    # Index symbols
    row = L.Symbol("row")
    col = L.Symbol("col")
    num = L.Symbol("num")

    combinations = L.Symbol("combinations")

    code = [
        L.Comment("Declare two dimensional array that holds combinations of derivatives and initialise"),
        L.ArrayDecl(int_type, combinations, (max_num_derivatives, max_degree)),
        L.ForRange(row, 0, max_num_derivatives, body=
            L.ForRange(col, 0, max_degree, body=
                L.Assign(combinations[row][col], 0))),
        L.Comment("Generate combinations of derivatives"),
        L.ForRange(row, 1, num_derivatives, body=
            L.ForRange(num, 0, row, body=
                L.For(L.VariableDecl(int_type, col, order-1),
                          L.GE(col, 0), L.PreDecrement(col), body=[
                    L.If(L.GT(combinations[row][col] + 1, tdim-1),
                        L.Assign(combinations[row][col], 0)),
                    L.Else([
                        L.AssignAdd(combinations[row][col], 1),
                        L.Break(),
                        ])
                    ])
                )
            ),
        ]

    # Python equivalent precomputed for each valid order:
    max_order = 3 # FIXME: degree of element?
    num_derivatives = max_num_derivatives
    #combinations = numpy.zeros((max_order, max_num_derivatives, max_degree))
    for order in range(1, max_order):
        combinations = numpy.zeros((num_derivatives, order))
        for row in range(1, num_derivatives):
            for num in range(0, row):
                for col in range(order-1, -1, -1):
                    if combinations[row][col] > tdim - 2:
                        combinations[row][col] = 0
                    else:
                        combinations[row][col] += 1
                        break

    return code


def generate_accumulation_code(L, dofs_data, num_derivatives, coefficients_for_dof):
    code = []
    # FIXME
    return code

def _compute_reference_derivatives(data, dof_data):  # FIXME: Convert this
    """Compute derivatives on the reference element by recursively multiply coefficients with
    the relevant derivatives of the polynomial base until the requested order of derivatives
    has been reached. After this take the dot product with the basisvalues."""

    # FIXME: Get symbols for these:
    f_basisvalues = format["basisvalues"]
    f_dmats = format["dmats"]
    f_dmats_old = format["dmats old"]
    f_derivatives = format["reference derivatives"]

    f_r, f_s, f_t, f_u = format["free indices"]

    tdim = data["topological_dimension"]
    gdim = data["geometric_dimension"]
    max_degree = data["max_degree"]

    # Get number of components.
    num_components = dof_data["num_components"]

    # Get shape of derivative matrix (they should all have the same shape) and
    # verify that it is a square matrix.
    shape_dmats = numpy.shape(dof_data["dmats"][0])
    if shape_dmats[0] != shape_dmats[1]:
        error("Something is wrong with the dmats:\n%s" % str(dof_data["dmats"]))

    code = [L.Comment("Compute reference derivatives.")]

    # Declare pointer to array that holds derivatives on the FIAT element
    code += [L.Comment("Declare array of derivatives on FIAT element.")]
    # The size of the array of reference derivatives is equal to the number of derivatives
    # times the number of components of the basis element
    if num_components == 1:
        num_vals = num_derivatives
    else:
        num_vals = num_components * num_derivatives

    nds = tdim**max_degree * num_components
    code += [L.ArrayDecl("double", f_derivatives, nds, values=0.0)]

    # Declare matrix of dmats (which will hold the matrix product of all combinations)
    # and dmats_old which is needed in order to perform the matrix product.
    eye = numpy.eye(shape_dmats[0])
    code += [
        L.Comment("Declare derivative matrix (of polynomial basis)."),
        L.ArrayDecl("double", f_dmats(""), shape_dmats, values=eye),
        L.Comment("Declare (auxiliary) derivative matrix (of polynomial basis)."),
        L.VariableDecl("double", f_dmats_old, shape_dmats, values=eye)
        ]

    # Compute dmats as a recursive matrix product
    # FIXME:
    dmats_lines = _compute_dmats(len(dof_data["dmats"]), shape_dmats, [f_s, f_t, f_u], f_r, _t)

    # Compute derivatives for all components
    lines_c = []
    for i in range(num_components):
        name = f_derivatives[i*num_derivatives + f_r]
        coeffs = coefficients_for_dof[i][f_s]
        dmats = f_dmats("")[f_s, f_t]
        basis = f_basisvalues[f_t]
        lines_c += [L.AssignAdd(name, coeffs * dmats * basis)]

    # Generate loop over number of derivatives.
    # Loop all derivatives and compute value of the derivative as:
    # deriv_on_ref[r] = coeff[dof][s]*dmat[s][t]*basis[t]
    code += [L.Comment("Loop possible derivatives.")]
    code += [L.ForRange(f_r, 0, num_derivatives, body=[
                L.StatementList(dmats_lines),
                L.ForRange(f_s, 0, shape_dmats[0], body=
                    L.ForRange(f_t, 0, shape_dmats[1], body=
                        lines_c)),
                ])
            ]
    return code




# dmats_old = eye
def _reset_dmats(shape_dmats, indices):
    "Set values in dmats equal to the identity matrix."
    f_assign = format["assign"]
    f_float = format["floating point"]
    i, j = indices

    code = [format["comment"]("Resetting dmats values to compute next derivative.")]
    dmats_old = format["component"](format["dmats"](""), [i, j])
    lines = [f_assign(dmats_old, f_float(0.0))]
    lines += [format["if"](i + format["is equal"] + j,
              f_assign(dmats_old, f_float(1.0)))]
    loop_vars = [(i, 0, shape_dmats[0]), (j, 0, shape_dmats[1])]
    code += format["generate loop"](lines, loop_vars)
    return code

# dmats_old = dmats
# dmats = 0.0
def _update_dmats(shape_dmats, indices):
    "Update values in dmats_old with values in dmats and set values in dmats to zero."
    f_assign = format["assign"]
    f_component = format["component"]
    i, j = indices

    code = [format["comment"]("Updating dmats_old with new values and resetting dmats.")]
    dmats = f_component(format["dmats"](""), [i, j])
    dmats_old = f_component(format["dmats old"], [i, j])
    lines = [f_assign(dmats_old, dmats),
             f_assign(dmats, format["floating point"](0.0))]
    loop_vars = [(i, 0, shape_dmats[0]), (j, 0, shape_dmats[1])]
    code += format["generate loop"](lines, loop_vars)
    return code


def _dmats_product(shape_dmats, index, i, indices):
    "Create product to update dmats."
    t, u = indices
    tu = t + u  # FIXME: Define index another way, this won't work with L.Symbol!

    # FIXME: Get the right symbols
    dmats = format["dmats"]("")
    dmats_i = format["dmats"](i)
    dmats_old = format["dmats old"]

    # dmats[t, u] = sum_k dmats_i[t, k] * dmats_old[k, u]
    theloop = [
        L.ForRange(t, 0, shape_dmats[0], body=
            L.ForRange(u, 0, shape_dmats[1], body=
                L.ForRange(tu, 0, shape_dmats[0], body=
                    L.AssignAdd(dmats[t, u], dmats_i[t, tu] * dmats_old[tu, u]))))
        ]

    code = [L.If(L.EQ(index, i), theloop)]
    return code


def _compute_dmats(num_dmats, shape_dmats, available_indices, deriv_index, _t):
    "Compute values of dmats as a matrix product."
    f_comment = format["comment"]
    s, t, u = available_indices

    # Reset dmats_old
# dmats_old = eye
    reset_dmats_code = _reset_dmats(shape_dmats, [t, u])

    # Set dmats matrix equal to dmats_old
# dmats_old = dmats
# dmats = 0.0
    update_dmats_code = _update_dmats(shape_dmats, [t, u])

    # Create dmats matrix by multiplication
    comb = format["component"](
        format["derivative combinations"](_t),
        [deriv_index, s]
        )
    products_code = []
    for i in range(num_dmats):
        products_code += _dmats_product(shape_dmats, comb, i, [t, u])

    code = [
        reset_dmats_code,
        L.Comment("Looping derivative order to generate dmats."),
        L.ForRange(s, 0, order, body=[
            update_dmats_code,
            L.Comment("Update dmats using an inner product."),
            products_code
            ])
        ]
    return code
