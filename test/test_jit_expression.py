# Copyright (C) 2019-2024 Michal Habera and Jørgen S. Dokken
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import basix
import basix.ufl
import cffi
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit


def test_matvec(compile_args):
    """Test evaluation of linear rank-0 form.

    Evaluates expression c * A_ij * f_j where c is a Constant,
    A_ij is a user specified constant matrix and f_j is j-th component
    of user specified vector-valued finite element function (in P1 space).

    """
    e = basix.ufl.element("P", "triangle", 1, shape=(2,))
    mesh = ufl.Mesh(e)
    V = ufl.FunctionSpace(mesh, e)
    f = ufl.Coefficient(V)

    a_mat = np.array([[1.0, 2.0], [1.0, 1.0]])
    a = ufl.as_matrix(a_mat)
    expr = ufl.Constant(mesh) * ufl.dot(a, f)

    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    obj, module, code = ffcx.codegeneration.jit.compile_expressions(
        [(expr, points)], cffi_extra_compile_args=compile_args
    )

    ffi = cffi.FFI()
    expression = obj[0]

    dtype = np.float64
    c_type = "double"
    xdtype = np.float64
    c_xtype = "double"

    A = np.zeros((3, 2), dtype=dtype)
    f_mat = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    # Coefficient storage XYXYXY
    w = np.array(f_mat.T.flatten(), dtype=dtype)
    c = np.array([0.5], dtype=dtype)
    entity_index = np.array([0], dtype=np.intc)
    quad_perm = np.array([0], dtype=np.dtype("uint8"))

    # Coords storage XYZXYZXYZ
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)
    expression.tabulate_tensor_float64(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", entity_index.ctypes.data),
        ffi.cast("uint8_t *", quad_perm.ctypes.data),
    )

    # Check the computation against correct NumPy value
    assert np.allclose(A, 0.5 * np.dot(a_mat, f_mat).T)

    # Prepare NumPy array of points attached to the expression
    length = expression.num_points * expression.entity_dimension
    points_kernel = np.frombuffer(
        ffi.buffer(expression.points, length * ffi.sizeof("double")), np.double
    )
    points_kernel = points_kernel.reshape(points.shape)
    assert np.allclose(points, points_kernel)

    # Check the value shape attached to the expression
    value_shape = np.frombuffer(
        ffi.buffer(expression.value_shape, expression.num_components * ffi.sizeof("int")), np.intc
    )
    assert np.allclose(expr.ufl_shape, value_shape)


def test_rank1(compile_args):
    """Tests evaluation of rank-1 form.

    Builds a linear operator which takes vector-valued functions in P1 space
    and evaluates expression [u_y, u_x] + grad(u_x) at specified points.

    """
    e = basix.ufl.element("P", "triangle", 1, shape=(2,))
    mesh = ufl.Mesh(e)

    V = ufl.FunctionSpace(mesh, e)
    u = ufl.TrialFunction(V)

    expr = ufl.as_vector([u[1], u[0]]) + ufl.grad(u[0])

    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    obj, module, code = ffcx.codegeneration.jit.compile_expressions(
        [(expr, points)], cffi_extra_compile_args=compile_args
    )

    ffi = cffi.FFI()
    expression = obj[0]

    dtype = np.float64
    c_type = "double"
    xdtype = np.float64
    c_xtype = "double"

    # 2 components for vector components of expression
    # 3 points of evaluation
    # 6 degrees of freedom for rank1 form
    A = np.zeros((3, 2, 6), dtype=dtype)

    # Coefficient storage XYXYXY
    w = np.array([0.0], dtype=dtype)
    c = np.array([0.0], dtype=dtype)
    entity_index = np.array([0], dtype=np.intc)
    quad_perm = np.array([0], dtype=np.dtype("uint8"))

    # Coords storage XYZXYZXYZ
    coords = np.zeros((points.shape[0], 3), dtype=xdtype)
    coords[:, :2] = points
    expression.tabulate_tensor_float64(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", entity_index.ctypes.data),
        ffi.cast("uint8_t *", quad_perm.ctypes.data),
    )

    f = np.array([[1.0, 2.0, 3.0], [-4.0, -5.0, 6.0]])

    # Apply the operator on some test input data
    u_ffcx = np.einsum("ijk,k", A, f.T.flatten())

    # Compute the correct values using NumPy
    # Gradf0 is gradient of f[0], each component of the gradient is constant
    gradf0 = np.array(
        [
            [f[0, 1] - f[0, 0], f[0, 1] - f[0, 0], f[0, 1] - f[0, 0]],
            [f[0, 2] - f[0, 0], f[0, 2] - f[0, 0], f[0, 2] - f[0, 0]],
        ]
    )

    u_correct = np.array([f[1], f[0]]) + gradf0

    assert np.allclose(u_ffcx, u_correct.T)


def test_elimiate_zero_tables_tensor(compile_args):
    """Test elimination of tensor-valued expressions with zero tables"""
    cell = "tetrahedron"
    c_el = basix.ufl.element("P", cell, 1, shape=(3,))
    mesh = ufl.Mesh(c_el)

    e = basix.ufl.element("P", cell, 1)
    V = ufl.FunctionSpace(mesh, e)
    u = ufl.Coefficient(V)
    expr = ufl.sym(ufl.as_tensor([[u, u.dx(0).dx(0), 0], [u.dx(1), u.dx(1), 0], [0, 0, 0]]))

    # Get vertices of cell
    # Coords storage XYZXYZXYZ
    basix_c_e = basix.create_element(
        basix.ElementFamily.P, basix.cell.string_to_type(cell), 1, discontinuous=False
    )
    coords = basix_c_e.points

    # Using same basix element for coordinate element and coefficient
    coeff_points = basix_c_e.points

    # Compile expression at interpolation points of second order Lagrange space
    b_el = basix.create_element(
        basix.ElementFamily.P, basix.cell.string_to_type(cell), 0, discontinuous=True
    )
    points = b_el.points
    obj, module, code = ffcx.codegeneration.jit.compile_expressions(
        [(expr, points)], cffi_extra_compile_args=compile_args
    )

    ffi = cffi.FFI()
    expression = obj[0]

    dtype = np.float64
    c_type = "double"
    c_xtype = "double"

    output = np.zeros(9 * points.shape[0], dtype=dtype)

    # Define coefficients for u = x + 2 * y
    u_coeffs = u_coeffs = coeff_points.T[0] + 2 * coeff_points.T[1]
    consts = np.array([], dtype=dtype)
    entity_index = np.array([0], dtype=np.intc)
    quad_perm = np.array([0], dtype=np.dtype("uint8"))

    expression.tabulate_tensor_float64(
        ffi.cast(f"{c_type} *", output.ctypes.data),
        ffi.cast(f"{c_type} *", u_coeffs.ctypes.data),
        ffi.cast(f"{c_type} *", consts.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", entity_index.ctypes.data),
        ffi.cast("uint8_t *", quad_perm.ctypes.data),
    )

    def exact_expr(x):
        val = np.zeros((9, x.shape[1]), dtype=dtype)
        val[0] = x[0] + 2 * x[1]
        val[1] = 0 + 0.5 * 2
        val[3] = 0.5 * 2 + 0
        val[4] = 2
        return val.T

    exact = exact_expr(points.T)

    assert np.allclose(exact, output)


def test_grad_constant(compile_args):
    """Test constant numbering.

    Test if numbering of constants are correct after UFL eliminates the
    constant inside the gradient.
    """
    c_el = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
    mesh = ufl.Mesh(c_el)

    x = ufl.SpatialCoordinate(mesh)
    first_constant = ufl.Constant(mesh)
    second_constant = ufl.Constant(mesh)
    expr = second_constant * ufl.Dx(x[0] ** 2 + first_constant, 0)

    dtype = np.float64
    points = np.array([[0.33, 0.25]], dtype=dtype)

    obj, _, _ = ffcx.codegeneration.jit.compile_expressions(
        [(expr, points)], cffi_extra_compile_args=compile_args
    )

    ffi = cffi.FFI()
    expression = obj[0]

    c_type = "double"
    c_xtype = "double"

    output = np.zeros(1, dtype=dtype)

    # Define constants
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=dtype)
    u_coeffs = np.array([], dtype=dtype)
    consts = np.array([3, 7], dtype=dtype)
    entity_index = np.array([0], dtype=np.intc)
    quad_perm = np.array([0], dtype=np.dtype("uint8"))

    expression.tabulate_tensor_float64(
        ffi.cast(f"{c_type} *", output.ctypes.data),
        ffi.cast(f"{c_type} *", u_coeffs.ctypes.data),
        ffi.cast(f"{c_type} *", consts.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", entity_index.ctypes.data),
        ffi.cast("uint8_t *", quad_perm.ctypes.data),
    )

    assert output[0] == pytest.approx(consts[1] * 2 * points[0, 0])


def test_facet_expression(compile_args):
    """Test facet expression containing a facet normal on a manifold."""
    c_el = basix.ufl.element("Lagrange", "triangle", 1, shape=(3,))
    mesh = ufl.Mesh(c_el)

    n = ufl.FacetNormal(mesh)
    expr = n

    dtype = np.float64
    points = np.array([[0.5]], dtype=dtype)

    obj, _, _ = ffcx.codegeneration.jit.compile_expressions(
        [(expr, points)], cffi_extra_compile_args=compile_args
    )

    ffi = cffi.FFI()
    expression = obj[0]

    c_type = "double"
    c_xtype = "double"

    output = np.zeros(3, dtype=dtype)

    # Define constants
    coords = np.array([[0.3, 0.6, 0.1], [1.2, 0.4, 0.2], [1.3, 1.4, 0.3]], dtype=dtype)
    u_coeffs = np.array([], dtype=dtype)
    consts = np.array([], dtype=dtype)
    entity_index = np.array([0], dtype=np.intc)
    quad_perm = np.array([0], dtype=np.dtype("uint8"))
    tangents = np.array([coords[1] - coords[2], coords[2] - coords[0], coords[0] - coords[1]])
    midpoints = np.array(
        [
            coords[1] + (coords[2] - coords[1]) / 2,
            coords[0] + (coords[2] - coords[0]) / 2,
            coords[1] + (coords[1] - coords[0]) / 2,
        ]
    )
    for i, (tangent, midpoint) in enumerate(zip(tangents, midpoints)):
        # normalize tangent
        tangent /= np.linalg.norm(tangent)
        # Tabulate facet normal
        output[:] = 0
        entity_index[0] = i
        expression.tabulate_tensor_float64(
            ffi.cast(f"{c_type} *", output.ctypes.data),
            ffi.cast(f"{c_type} *", u_coeffs.ctypes.data),
            ffi.cast(f"{c_type} *", consts.ctypes.data),
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.cast("int *", entity_index.ctypes.data),
            ffi.cast("uint8_t *", quad_perm.ctypes.data),
        )
        # Assert that facet normal is perpendicular to tangent
        assert np.isclose(np.dot(output, tangent), 0)

        # Check that norm of facet normal is 1
        assert np.isclose(np.linalg.norm(output), 1)

        # Check that facet normal is pointing out of the cell
        assert np.dot(midpoint - coords[i], output) > 0
