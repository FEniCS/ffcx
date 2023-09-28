# -*- coding: utf-8 -*-
# Copyright (C) 2019-2022 Michal Habera and Jørgen S. Dokken
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import cffi
import numpy as np

import basix
import basix.ufl
import ffcx.codegeneration.jit
import ufl
from ffcx.codegeneration.utils import cdtype_to_numpy, scalar_to_value_type


def float_to_type(name):
    """Map a string name to C and NumPy types"""
    if name == "double":
        return "double", np.float64
    elif name == "double complex":
        return "double _Complex", np.complex128
    elif name == "float":
        return "float", np.float32
    elif name == "float complex":
        return "float _Complex", np.complex64
    elif name == "long double":
        return "long double", np.longdouble
    else:
        raise RuntimeError("Unknown C type for: {}".format(name))


def test_matvec(compile_args):
    """Test evaluation of linear rank-0 form.

    Evaluates expression c * A_ij * f_j where c is a Constant,
    A_ij is a user specified constant matrix and f_j is j-th component
    of user specified vector-valued finite element function (in P1 space).

    """
    e = basix.ufl.element("P", "triangle", 1, shape=(2, ))
    mesh = ufl.Mesh(e)
    V = ufl.FunctionSpace(mesh, e)
    f = ufl.Coefficient(V)

    a_mat = np.array([[1.0, 2.0], [1.0, 1.0]])
    a = ufl.as_matrix(a_mat)
    expr = ufl.Constant(mesh) * ufl.dot(a, f)

    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    obj, module, code = ffcx.codegeneration.jit.compile_expressions(
        [(expr, points)], cffi_extra_compile_args=compile_args)

    ffi = cffi.FFI()
    expression = obj[0]

    c_type, np_type = float_to_type("double")

    A = np.zeros((3, 2), dtype=np_type)
    f_mat = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    # Coefficient storage XYXYXY
    w = np.array(f_mat.T.flatten(), dtype=np_type)
    c = np.array([0.5], dtype=np_type)
    entity_index = np.array([0], dtype=np.intc)
    quad_perm = np.array([0], dtype=np.dtype("uint8"))

    geom_type = scalar_to_value_type(c_type)
    np_gtype = cdtype_to_numpy(geom_type)

    # Coords storage XYZXYZXYZ
    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np_gtype)
    expression.tabulate_tensor_float64(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast(f'{geom_type} *', coords.ctypes.data),
        ffi.cast('int *', entity_index.ctypes.data),
        ffi.cast('uint8_t *', quad_perm.ctypes.data))

    # Check the computation against correct NumPy value
    assert np.allclose(A, 0.5 * np.dot(a_mat, f_mat).T)

    # Prepare NumPy array of points attached to the expression
    length = expression.num_points * expression.topological_dimension
    points_kernel = np.frombuffer(ffi.buffer(expression.points, length * ffi.sizeof("double")), np.double)
    points_kernel = points_kernel.reshape(points.shape)
    assert np.allclose(points, points_kernel)

    # Check the value shape attached to the expression
    value_shape = np.frombuffer(ffi.buffer(expression.value_shape,
                                           expression.num_components * ffi.sizeof("int")), np.intc)
    assert np.allclose(expr.ufl_shape, value_shape)


def test_rank1(compile_args):
    """Tests evaluation of rank-1 form.

    Builds a linear operator which takes vector-valued functions in P1 space
    and evaluates expression [u_y, u_x] + grad(u_x) at specified points.

    """
    e = basix.ufl.element("P", "triangle", 1, shape=(2, ))
    mesh = ufl.Mesh(e)

    V = ufl.FunctionSpace(mesh, e)
    u = ufl.TrialFunction(V)

    expr = ufl.as_vector([u[1], u[0]]) + ufl.grad(u[0])

    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    obj, module, code = ffcx.codegeneration.jit.compile_expressions(
        [(expr, points)], cffi_extra_compile_args=compile_args)

    ffi = cffi.FFI()
    expression = obj[0]

    c_type, np_type = float_to_type("double")
    geom_type = scalar_to_value_type(c_type)
    np_gtype = cdtype_to_numpy(geom_type)

    # 2 components for vector components of expression
    # 3 points of evaluation
    # 6 degrees of freedom for rank1 form
    A = np.zeros((3, 2, 6), dtype=np_type)

    # Coefficient storage XYXYXY
    w = np.array([0.0], dtype=np_type)
    c = np.array([0.0], dtype=np_type)
    entity_index = np.array([0], dtype=np.intc)
    quad_perm = np.array([0], dtype=np.dtype("uint8"))

    # Coords storage XYZXYZXYZ
    coords = np.zeros((points.shape[0], 3), dtype=np_gtype)
    coords[:, :2] = points
    expression.tabulate_tensor_float64(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast(f'{geom_type} *', coords.ctypes.data),
        ffi.cast('int *', entity_index.ctypes.data),
        ffi.cast('uint8_t *', quad_perm.ctypes.data))

    f = np.array([[1.0, 2.0, 3.0], [-4.0, -5.0, 6.0]])

    # Apply the operator on some test input data
    u_ffcx = np.einsum("ijk,k", A, f.T.flatten())

    # Compute the correct values using NumPy
    # Gradf0 is gradient of f[0], each component of the gradient is constant
    gradf0 = np.array([[f[0, 1] - f[0, 0], f[0, 1] - f[0, 0], f[0, 1] - f[0, 0]],
                       [f[0, 2] - f[0, 0], f[0, 2] - f[0, 0], f[0, 2] - f[0, 0]]])

    u_correct = np.array([f[1], f[0]]) + gradf0

    assert np.allclose(u_ffcx, u_correct.T)


def test_elimiate_zero_tables_tensor(compile_args):
    """
    Test elimination of tensor-valued expressions with zero tables
    """
    cell = "tetrahedron"
    c_el = basix.ufl.element("P", cell, 1, shape=(3, ))
    mesh = ufl.Mesh(c_el)

    e = basix.ufl.element("P", cell, 1)
    V = ufl.FunctionSpace(mesh, e)
    u = ufl.Coefficient(V)
    expr = ufl.sym(ufl.as_tensor([[u, u.dx(0).dx(0), 0],
                                  [u.dx(1), u.dx(1), 0],
                                  [0, 0, 0]]))

    # Get vertices of cell
    # Coords storage XYZXYZXYZ
    basix_c_e = basix.create_element(basix.ElementFamily.P,
                                     basix.cell.string_to_type(cell), 1,
                                     discontinuous=False)
    coords = basix_c_e.points

    # Using same basix element for coordinate element and coefficient
    coeff_points = basix_c_e.points

    # Compile expression at interpolation points of second order Lagrange space
    b_el = basix.create_element(basix.ElementFamily.P, basix.cell.string_to_type(cell),
                                0, discontinuous=True)
    points = b_el.points
    obj, module, code = ffcx.codegeneration.jit.compile_expressions(
        [(expr, points)], cffi_extra_compile_args=compile_args)

    ffi = cffi.FFI()
    expression = obj[0]
    c_type, np_type = float_to_type("double")

    output = np.zeros(9 * points.shape[0], dtype=np_type)

    # Define coefficients for u = x + 2 * y
    u_coeffs = u_coeffs = coeff_points.T[0] + 2 * coeff_points.T[1]
    consts = np.array([], dtype=np_type)
    entity_index = np.array([0], dtype=np.intc)
    quad_perm = np.array([0], dtype=np.dtype("uint8"))

    expression.tabulate_tensor_float64(
        ffi.cast('{type} *'.format(type=c_type), output.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), u_coeffs.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), consts.ctypes.data),
        ffi.cast('double *', coords.ctypes.data),
        ffi.cast('int *', entity_index.ctypes.data),
        ffi.cast('uint8_t *', quad_perm.ctypes.data))

    def exact_expr(x):
        val = np.zeros((9, x.shape[1]), dtype=np_type)
        val[0] = x[0] + 2 * x[1]
        val[1] = 0 + 0.5 * 2
        val[3] = 0.5 * 2 + 0
        val[4] = 2
        return val.T
    exact = exact_expr(points.T)

    assert np.allclose(exact, output)
