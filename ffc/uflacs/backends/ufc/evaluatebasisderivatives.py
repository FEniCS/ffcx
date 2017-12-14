# -*- coding: utf-8 -*-
"""Work in progress translation of FFC evaluatebasisderivatives code to uflacs CNodes format."""

import numpy
import math

from ffc.log import error
from ffc.uflacs.backends.ufc.utils import generate_error
from ffc.uflacs.backends.ufc.evaluatebasis import _generate_compute_basisvalues, tabulate_coefficients
from ffc.uflacs.backends.ufc.evalderivs import _generate_combinations
from ffc.uflacs.backends.ufc.jacobian import jacobian, inverse_jacobian, orientation, fiat_coordinate_mapping, _mapping_transform

# Used for various indices and arrays in this file
index_type = "std::size_t"

def _compute_reference_derivatives(L, data, dof_data):
    """Compute derivatives on the reference element by recursively
    multiply coefficients with the relevant derivatives of the
    polynomial base until the requested order of derivatives has been
    reached. After this take the dot product with the basisvalues.

    """

    # Prefetch formats to speed up code generation

    tdim = data["topological_dimension"]
    gdim = data["geometric_dimension"]
    max_degree = data["max_degree"]

    _t, _g = ("","") if (tdim == gdim) else ("_t", "_g")
    num_derivs_t = L.Symbol("num_derivatives" + _t)
    num_derivs_g = L.Symbol("num_derivatives" + _g)

    # Get number of components.
    num_components = dof_data["num_components"]

    # Get shape of derivative matrix (they should all have the same
    # shape) and verify that it is a square matrix.
    shape_dmats = dof_data["dmats"][0].shape
    if shape_dmats[0] != shape_dmats[1]:
        error("Something is wrong with the dmats:\n%s" % str(dof_data["dmats"]))

    code = [L.Comment("Compute reference derivatives.")]

    # Declare pointer to array that holds derivatives on the FIAT element
    code += [L.Comment("Declare array of derivatives on FIAT element.")]
    # The size of the array of reference derivatives is equal to the
    # number of derivatives times the number of components of the
    # basis element
    num_vals = num_components*num_derivs_t
    nds = tdim**max_degree * num_components

    mapping = dof_data["mapping"]
    if "piola" in mapping and "double" not in mapping and gdim > num_components :
        # In either of the Piola cases, the value space of the
        # derivatives is the geometric dimension rather than the
        # topological dimension.  Increase size of derivatives array
        # if needed.
        nds = tdim**max_degree * gdim

    derivatives = L.Symbol("derivatives")
    code += [L.ArrayDecl("double", derivatives, nds, 0.0)]

    # Declare matrix of dmats (which will hold the matrix product of
    # all combinations) and dmats_old which is needed in order to
    # perform the matrix product.
    code += [L.Comment("Declare derivative matrix (of polynomial basis).")]
    value = numpy.eye(shape_dmats[0])
    dmats = L.Symbol("dmats")
    code += [L.ArrayDecl("double", dmats, shape_dmats, value)]
    code += [L.Comment("Declare (auxiliary) derivative matrix (of polynomial basis).")]
    dmats_old = L.Symbol("dmats_old")
    code += [L.ArrayDecl("double", dmats_old, shape_dmats, value)]

    r = L.Symbol("r")
    s = L.Symbol("s")
    t = L.Symbol("t")
    n = L.Symbol("n")

    lines = [L.Comment("Reset dmats to identity"),
             L.MemZero(L.AddressOf(dmats[0][0]), shape_dmats[0]*shape_dmats[1]),
             L.ForRange(t, 0, shape_dmats[0], index_type=index_type, body=[L.Assign(dmats[t][t], 1.0)])]

    lines_s = [L.MemCopy(L.AddressOf(dmats[0][0]), L.AddressOf(dmats_old[0][0]), shape_dmats[0]*shape_dmats[1]),
               L.MemZero(L.AddressOf(dmats[0][0]), shape_dmats[0]*shape_dmats[1])]

    lines_s += [L.Comment("Update dmats using an inner product.")]

    # Create dmats matrix by multiplication
    comb = L.Symbol("combinations" + _t)

    for i in range(len(dof_data["dmats"])):
        lines_s += [L.Comment("_dmats_product(shape_dmats, comb[r][s], %d)" % i)]
        dmats_i = L.Symbol("dmats%d" % i)
        u = L.Symbol("u")
        tu = L.Symbol("tu")
        lines_comb = [L.ForRange(t, 0, shape_dmats[0], index_type=index_type,
              body = [L.ForRange(u, 0, shape_dmats[1], index_type=index_type,
              body = [L.ForRange(tu, 0, shape_dmats[1], index_type=index_type,
              body = [L.AssignAdd(dmats[t][u], dmats_old[tu][u] * dmats_i[t][tu])])])])]

        lines_s += [L.If(L.EQ(comb[n - 1][r][s], i), lines_comb)]

    lines += [L.Comment("Looping derivative order to generate dmats."),
              L.ForRange(s, 0, n, index_type=index_type, body=lines_s)]

    # Compute derivatives for all components
    lines_c = []
    basisvalues = L.Symbol("basisvalues")
    for i in range(num_components):
        coeffs = L.Symbol("coefficients%d" % i)
        lines_c += [L.AssignAdd(derivatives[i*num_derivs_t + r], coeffs[s]*dmats[s,t]*basisvalues[t])]
    lines += [L.ForRange(s, 0, shape_dmats[0], index_type=index_type,
        body=[L.ForRange(t, 0, shape_dmats[1], index_type=index_type, body=lines_c)])]

    lines += _mapping_transform(L, data, dof_data, derivatives, r, num_derivs_t)

    code += [L.Comment("Loop possible derivatives."),
             L.ForRange(r, 0, num_derivs_t, index_type=index_type, body=lines)]

    return code


def _transform_derivatives(L, data, dof_data):
    """Transform derivatives back to the physical element by applying the
    transformation matrix."""

    tdim = data["topological_dimension"]
    gdim = data["geometric_dimension"]

    _t, _g = ("","") if (tdim == gdim) else ("_t", "_g")
    num_derivs_t = L.Symbol("num_derivatives" + _t)
    num_derivs_g = L.Symbol("num_derivatives" + _g)

    # Get number of components and offset.
    num_components = dof_data["num_components"]
    reference_offset = dof_data["reference_offset"]
    physical_offset = dof_data["physical_offset"]
    offset = reference_offset  # physical_offset # FIXME: Should be physical offset but that breaks tests

    mapping = dof_data["mapping"]
    if "piola" in mapping and "double" not in mapping:
        # In either of the Piola cases, the value space of the derivatives
        # is the geometric dimension rather than the topological dimension.
        num_components = gdim

    code = [L.Comment("Transform derivatives back to physical element")]

    lines = []
    r = L.Symbol("r")
    s = L.Symbol("s")
    values = L.Symbol("values")
    transform = L.Symbol("transform")
    derivatives = L.Symbol("derivatives")
    for i in range(num_components):
        lines += [L.AssignAdd(values[(offset + i)*num_derivs_g + r],
                              transform[r, s]*derivatives[i*num_derivs_t + s])]

    code += [L.ForRange(r, 0, num_derivs_g, index_type=index_type,
                        body=[L.ForRange(s, 0, num_derivs_t, index_type=index_type, body=lines)])]

    return code


def _tabulate_dmats(L, dof_data):
    "Tabulate the derivatives of the polynomial base"

    # Get derivative matrices (coefficients) of basis functions,
    # computed by FIAT at compile time.

    code = [L.Comment("Tables of derivatives of the polynomial base (transpose).")]

    # Generate tables for each spatial direction.
    for i, dmat in enumerate(dof_data["dmats"]):

        # Extract derivatives for current direction (take transpose,
        # FIAT_NEW PolynomialSet.tabulate()).
        matrix = numpy.transpose(dmat)

        # Get shape and check dimension (This is probably not needed).
        #        shape = numpy.shape(matrix)
        if not (matrix.shape[0] == matrix.shape[1] == dof_data["num_expansion_members"]):
            error("Something is wrong with the shape of dmats.")

        # Declare varable name for coefficients.
        table = L.Symbol("dmats%d" % i)
        code += [L.ArrayDecl("static const double", table, matrix.shape, matrix)]

    return code


def _generate_dof_code(L, data, dof_data):
    "Generate code for a basis."

    basisvalues = L.Symbol("basisvalues")
    Y = L.Symbol("Y")
    element_cellname = data["cellname"]
    embedded_degree = dof_data["embedded_degree"]
    num_members = dof_data["num_expansion_members"]
    code = _generate_compute_basisvalues(L, basisvalues, Y, element_cellname, embedded_degree, num_members)

    # Tabulate coefficients.
    code += tabulate_coefficients(L, dof_data)

    # Tabulate coefficients for derivatives.
    code += _tabulate_dmats(L, dof_data)

    # Compute the derivatives of the basisfunctions on the reference
    # (FIAT) element, as the dot product of the new coefficients and
    # basisvalues.
    code += _compute_reference_derivatives(L, data, dof_data)

    # Transform derivatives to physical element by multiplication with
    # the transformation matrix.
    code += _transform_derivatives(L, data, dof_data)

    return code


def _generate_transform(L, element_cellname, gdim, tdim, max_degree):
    """Generate the transformation matrix, which is used to transform
    derivatives from reference element back to the physical
    element.

    """

    max_g_d = gdim**max_degree
    max_t_d = tdim**max_degree

    K = L.Symbol("K")
    transform = L.Symbol("transform")
    col = L.Symbol("col")
    row = L.Symbol("row")

    _t, _g = ("","") if (tdim == gdim) else ("_t", "_g")
    comb_t = L.Symbol("combinations" + _t)
    comb_g = L.Symbol("combinations" + _g)
    num_derivatives_t = L.Symbol("num_derivatives" + _t)
    num_derivatives_g = L.Symbol("num_derivatives" + _g)

    K = L.FlattenedArray(K, dims=(tdim, gdim))

    code = [L.Comment("Declare transformation matrix")]
    code += [L.ArrayDecl("double", transform, (max_g_d, max_t_d), numpy.ones((max_g_d, max_t_d)))]
    code += [L.Comment("Construct transformation matrix")]
    k = L.Symbol("k")
    n = L.Symbol("n")
    inner_loop = L.ForRange(col, 0, num_derivatives_t, index_type=index_type,
           body=[L.ForRange(k, 0, n, index_type=index_type,
           body=[L.AssignMul(transform[row][col], K[comb_t[n-1][col][k], comb_g[n-1][row][k]])])])

    code += [L.ForRange(row, 0, num_derivatives_g, index_type=index_type, body=inner_loop)]

    return code


def generate_evaluate_basis_derivatives(L, data):
    """Evaluate the derivatives of an element basisfunction at a
    point. The values are computed as in FIAT as the matrix product of
    the coefficients (computed at compile time), basisvalues which are
    dependent on the coordinate and thus have to be computed at run
    time and combinations (depending on the order of derivative) of
    dmats tables which hold the derivatives of the expansion
    coefficients.

    """

    if isinstance(data, str):
        msg = "evaluate_basis_derivatives: %s" % data
        return [generate_error(L, msg, True)]

    # Initialise return code.
    code = []

    # Get the element cell domain, geometric and topological dimension.
    element_cellname = data["cellname"]
    gdim = data["geometric_dimension"]
    tdim = data["topological_dimension"]
    max_degree = data["max_degree"]
    physical_value_size = data["physical_value_size"]

    # Compute number of derivatives that has to be computed, and
    # declare an array to hold the values of the derivatives on the
    # reference element.
    values = L.Symbol("values")
    n = L.Symbol("n")
    dofs = L.Symbol("dofs")
    x = L.Symbol("x")
    coordinate_dofs = L.Symbol("coordinate_dofs")
    cell_orientation = L.Symbol("cell_orientation")

    if tdim == gdim:
        _t, _g = ("", "")
        num_derivatives_t = L.Symbol("num_derivatives")
        num_derivatives_g = L.Symbol("num_derivatives")
        code += [L.VariableDecl("std::size_t", num_derivatives_t, L.Call("std::pow", (tdim, n)))]
    else:
        _t, _g = ("_t", "_g")
        num_derivatives_t = L.Symbol("num_derivatives_t")
        num_derivatives_g = L.Symbol("num_derivatives_g")
        if max_degree > 0:
            code += [L.VariableDecl("std::size_t", num_derivatives_t, L.Call("std::pow", (tdim, n)))]
        code += [L.VariableDecl("std::size_t", num_derivatives_g, L.Call("std::pow", (gdim, n)))]

    # Reset all values.
    code += [L.MemZero(values, physical_value_size*num_derivatives_g)]

    # Handle values of argument 'n'.
    code += [L.Comment("Call evaluate_basis_all if order of derivatives is equal to zero.")]
    code += [L.If(L.EQ(n, 0), [L.Call("evaluate_basis", (L.Symbol("i"), values, x, coordinate_dofs, cell_orientation)), L.Return()])]

    # If max_degree is zero, return code (to avoid declarations such as
    # combinations[1][0]) and because there's nothing to compute.)
    if max_degree == 0:
        return code

    code += [L.Comment("If order of derivatives is greater than the maximum polynomial degree, return zeros.")]
    code += [L.If(L.GT(n, max_degree), [L.Return()])]

    # Generate geo code.
    code += jacobian(L, gdim, tdim, element_cellname)
    code += inverse_jacobian(L, gdim, tdim, element_cellname)
    if data["needs_oriented"] and tdim != gdim:
        code += orientation(L)

    code += fiat_coordinate_mapping(L, element_cellname, gdim)

    # Generate all possible combinations of derivatives.
    combinations_code_t, combinations_t = _generate_combinations(L, tdim, max_degree, n, num_derivatives_t, _t)
    code += combinations_code_t
    if tdim != gdim:
        combinations_code_g, combinations_g = _generate_combinations(L, gdim, max_degree, n, num_derivatives_g, "_g")
        code += combinations_code_g

    # Generate the transformation matrix.
    code += _generate_transform(L, element_cellname, gdim, tdim, max_degree)

    # Create code for all basis values (dofs).
    dof_cases = []
    for i, dof_data in enumerate(data["dofs_data"]):
        dof_cases.append((i, _generate_dof_code(L, data, dof_data)))
    code += [L.Switch(L.Symbol("i"), dof_cases)]
    return code

def generate_evaluate_basis_derivatives_all(L, data):
    """Like evaluate_basis, but return the values of all basis
    functions (dofs)."""

    if isinstance(data, str):
        msg = "evaluate_basis_derivatives_all: %s" % data
        return [generate_error(L, msg, True)]

    # Initialise return code
    code = []

    # FIXME: KBO: Figure out which return format to use, either:
    # [dN0[0]/dx, dN0[0]/dy, dN0[1]/dx, dN0[1]/dy, dN1[0]/dx,
    # dN1[0]/dy, dN1[1]/dx, dN1[1]/dy, ...]
    # or
    # [dN0[0]/dx, dN1[0]/dx, ..., dN0[1]/dx, dN1[1]/dx, ...,
    # dN0[0]/dy, dN1[0]/dy, ..., dN0[1]/dy, dN1[1]/dy, ...]
    # or
    # [dN0[0]/dx, dN0[1]/dx, ..., dN1[0]/dx, dN1[1]/dx, ...,
    # dN0[0]/dy, dN0[1]/dy, ..., dN1[0]/dy, dN1[1]/dy, ...]
    # for vector (tensor elements), currently returning option 1.

    # FIXME: KBO: For now, just call evaluate_basis_derivatives and
    # map values accordingly, this will keep the amount of code at a
    # minimum. If it turns out that speed is an issue (overhead from
    # calling evaluate_basis), we can easily generate all the code.

    # Get total value shape and space dimension for entire element
    # (possibly mixed).
    physical_value_size = data["physical_value_size"]
    space_dimension = data["space_dimension"]
    max_degree = data["max_degree"]
    gdim = data["geometric_dimension"]
    tdim = data["topological_dimension"]

    n = L.Symbol("n")
    x = L.Symbol("x")
    coordinate_dofs = L.Symbol("coordinate_dofs")
    cell_orientation = L.Symbol("cell_orientation")
    values = L.Symbol("values")

    # Special case where space dimension is one (constant elements).
    if space_dimension == 1:
        code += [L.Comment("Element is constant, calling evaluate_basis_derivatives.")]
        code += [L.Call("evaluate_basis_derivatives", (0, n, values, x, coordinate_dofs, cell_orientation))]
        return code

    # Compute number of derivatives.
    if tdim == gdim:
        num_derivatives = L.Symbol("num_derivatives")
    else:
        num_derivatives = L.Symbol("num_derivatives_g")

    # If n == 0, call evaluate_basis.
    code += [L.Comment("Call evaluate_basis_all if order of derivatives is equal to zero.")]
    code += [L.If(L.EQ(n, 0), [L.Call("evaluate_basis_all", (values, x, coordinate_dofs, cell_orientation)), L.Return()])]

    code += [L.VariableDecl("unsigned int", num_derivatives, L.Call("std::pow", (gdim, n)))]

    num_vals = physical_value_size * num_derivatives

    # Reset values.
    code += [L.Comment("Set values equal to zero.")]
    code += [L.MemZero(values, num_vals*space_dimension)]

    # If n > max_degree, return zeros.
    code += [L.Comment("If order of derivatives is greater than the maximum polynomial degree, return zeros.")]
    code += [L.If(L.GT(n, max_degree), [L.Return()])]

    # Declare helper value to hold single dof values and reset.
    code += [L.Comment("Helper variable to hold values of a single dof.")]
    nds = gdim**max_degree * physical_value_size
    dof_values = L.Symbol("dof_values")
    code += [L.ArrayDecl("double", dof_values, (nds,), 0.0)]

    # Create loop over dofs that calls evaluate_basis_derivatives for
    # a single dof and inserts the values into the global array.
    code += [L.Comment("Loop dofs and call evaluate_basis_derivatives.")]

    values = L.FlattenedArray(values, dims=(space_dimension, num_vals))
    r = L.Symbol("r")
    s = L.Symbol("s")
    loop_s = L.ForRange(s, 0, num_vals, index_type=index_type, body=[L.Assign(values[r, s], dof_values[s])])

    code += [L.ForRange(r, 0, space_dimension, index_type=index_type, body=[L.Call("evaluate_basis_derivatives", (r, n, dof_values, x, coordinate_dofs, cell_orientation)),
                                                                            loop_s])]
    return code
