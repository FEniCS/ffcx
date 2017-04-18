# -*- coding: utf-8 -*-
"""Work in progress translation of FFC evaluatebasisderivatives code to uflacs CNodes format."""

from six import string_types
import numpy
import math

from ffc.log import error
from ffc.uflacs.backends.ufc.utils import generate_error

# Used for various indices and arrays in this file
index_type = "std::size_t"

def _x_evaluate_basis_derivatives_all(L, data):
    """Like evaluate_basis, but return the values of all basis
    functions (dofs)."""

    if isinstance(data, string_types):
        return format["exception"]("evaluate_basis_derivatives_all: %s" % data)

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
        _g = ""
    else:
        _g = "_g"

    num_derivatives = L.Symbol("num_derivatives"+_g)

    # If n == 0, call evaluate_basis.
    code += [L.Comment("Call evaluate_basis_all if order of derivatives is equal to zero.")]
    code += [L.If(L.EQ(n, 0), [L.Call("evaluate_basis_all", (values, x, coordinate_dofs, cell_orientation)), L.Return()])]

    code += [L.Assign(num_derivatives, L.Call("std::pow", (gdim, n)))]

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

    # Create loop over dofs that calls evaluate_basis_derivatives for a single dof and
    # inserts the values into the global array.
    code += [L.Comment("Loop dofs and call evaluate_basis_derivatives.")]

    values = L.FlattenedArray(values, dims=(space_dimension, num_vals))
    r = L.Symbol("r")
    s = L.Symbol("s")
    loop_s = L.ForRange(s, 0, num_vals, index_type=index_type, body=[L.Assign(values[r, s], dof_values[s])])

    code += [L.ForRange(r, 0, space_dimension, index_type=index_type, body=[L.Call("evaluate_basis_derivatives", (r, n, dof_values, x, coordinate_dofs, cell_orientation)),
                                                                            loop_s])]
    return code
