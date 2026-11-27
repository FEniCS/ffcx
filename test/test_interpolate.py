# Copyright (C) 2026 Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
import itertools
import sys

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.analysis
import ffcx.codegeneration.jit
import ffcx.ir.representation
import ffcx.options
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("element", [("N1curl", {}), ("Lagrange", {"shape": (2,)})])
def test_interpolate(compile_args, dtype, element):

    cell = "triangle"
    family, el_kwargs = element

    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element = basix.ufl.element(family, cell, 2, **el_kwargs)
    V_int = ufl.FunctionSpace(domain, element)

    # Space containing other coefficients
    Q = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", cell, 2))
    w = ufl.Coefficient(Q)
    z = ufl.Coefficient(Q)
    q = ufl.Coefficient(Q)

    x = ufl.SpatialCoordinate(domain)
    c = ufl.Constant(domain)
    f = ufl.as_vector((z, -x[0])) + ufl.as_vector((q, q))
    If = ufl.Interpolate(f, V_int)
    J = (w * If[0] + If[1]) * ufl.dx
    f_ref = ufl.as_vector((x[1], -x[0])) + ufl.as_vector((1, 1))
    J_ref = (w * f_ref[0] + f_ref[1]) * ufl.dx

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    scale = 4.2
    coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, scale, 0.0]], dtype=xdtype).flatten()
    c = np.array([2.3], dtype=dtype)
    # Coefficients are ordered according to when they were created, thus
    # w, z, q
    # We set z to x[1], w, q to 1 and w to a some non-zero value
    q_size = Q.ufl_element().basix_element.dim
    assert q_size == 6
    d = np.empty(3 * q_size, dtype=dtype)
    d[0:q_size] = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]  # w
    d[q_size : 2 * q_size] = [0.0, 0.0, scale * 1.0, 0.5 * scale, 0.5 * scale, 0.0]  # z
    d[2 * q_size : 3 * q_size] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # q

    # Get kernel
    ffi = module.ffi
    form = compiled_forms[0]
    offsets = form.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1

    default_integral = form.form_integrals[offsets[cell]]

    A = np.zeros(1, dtype=dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", d.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    d_ref = np.zeros(q_size, dtype=dtype)
    d_ref[:] = d[0:q_size]  # w
    c_ref = c.copy()
    form_ref = compiled_forms[1]
    offsets_ref = form_ref.form_integral_offsets
    assert offsets_ref[cell + 1] - offsets_ref[cell] == 1
    ref_integral = form_ref.form_integrals[offsets_ref[cell]]

    A_ref = np.zeros(1, dtype=dtype)
    ref_kernel = getattr(ref_integral, f"tabulate_tensor_{dtype}")
    ref_kernel(
        ffi.cast(f"{c_type} *", A_ref.ctypes.data),
        ffi.cast(f"{c_type} *", d_ref.ctypes.data),
        ffi.cast(f"{c_type} *", c_ref.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    tol = np.finfo(dtype).eps * 100
    np.testing.assert_allclose(A, A_ref, atol=tol)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize(
    "elements",
    [
        (("Lagrange", 2, {}), ("Lagrange", 1, {})),
        (("Lagrange", 2, {"shape": (2,)}), ("Lagrange", 1, {"shape": (2,)})),
        (("N1curl", 2, {}), ("N1curl", 1, {})),
    ],
)
def test_interpolate_argument(compile_args, dtype, elements):
    """Interpolating an argument tabulates the target space against its dofs.

    ``I(phi_j) = sum_i M_ij psi_i`` for the reference interpolation operator
    ``M``, so the interpolated tensor must be ``M.T`` applied to the tensor of
    the same form posed on the target space.
    """
    cell = "triangle"
    (fam_a, deg_a, kw_a), (fam_t, deg_t, kw_t) = elements

    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element_a = basix.ufl.element(fam_a, cell, deg_a, **kw_a)
    element_t = basix.ufl.element(fam_t, cell, deg_t, **kw_t)
    V_a = ufl.FunctionSpace(domain, element_a)
    V_t = ufl.FunctionSpace(domain, element_t)

    def sum_components(f):
        """Contract all value components, so scalar and vector spaces are alike."""
        return f if f.ufl_shape == () else sum(f[i] for i in range(f.ufl_shape[0]))

    v_a, v_t = ufl.TrialFunction(V_a), ufl.TrialFunction(V_t)
    J = sum_components(ufl.Interpolate(v_a, V_t)) * ufl.dx
    J_ref = sum_components(v_t) * ufl.dx

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    cell_type = module.lib.cell

    def tabulate(form, size):
        offsets = form.form_integral_offsets
        assert offsets[cell_type + 1] - offsets[cell_type] == 1
        integral = form.form_integrals[offsets[cell_type]]
        A = np.zeros(size, dtype=dtype)
        kernel = getattr(integral, f"tabulate_tensor_{dtype}")
        kernel(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    A = tabulate(compiled_forms[0], element_a.dim)
    A_ref = tabulate(compiled_forms[1], element_t.dim)

    # Blocked elements interpolate each block component independently.
    def base_element(element):
        return element.sub_elements[0] if element.block_size > 1 else element

    interpolation = basix.compute_interpolation_operator(
        base_element(element_a).basix_element, base_element(element_t).basix_element
    )
    interpolation = np.kron(interpolation, np.eye(element_a.block_size))

    tol = np.finfo(dtype).eps * 100
    assert np.abs(A).max() > tol
    np.testing.assert_allclose(A, interpolation.T @ A_ref, atol=tol)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("interpolated", ["test", "trial", "both"])
@pytest.mark.parametrize(
    "family, kwargs", [("Lagrange", {}), ("Lagrange", {"shape": (2,)}), ("N1curl", {})]
)
def test_interpolate_argument_bilinear(compile_args, dtype, interpolated, family, kwargs):
    """Interpolating the arguments of a bilinear form.

    Each interpolated argument maps its own axis of the tensor through the
    reference interpolation operator, ``A = M_0^T A_ref M_1`` for the tensor
    ``A_ref`` of the same form posed on the target space, with ``M_i`` the
    identity for an argument that is not interpolated. The test and trial
    functions come from different spaces, so interpolating both gives two
    different tables into the one target space.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element_0 = basix.ufl.element(family, cell, 3, **kwargs)
    element_1 = basix.ufl.element(family, cell, 2, **kwargs)
    element_t = basix.ufl.element(family, cell, 1, **kwargs)
    V_0, V_1, V_t = (ufl.FunctionSpace(domain, e) for e in (element_0, element_1, element_t))

    v, u = ufl.TestFunction(V_0), ufl.TrialFunction(V_1)
    v_t, u_t = ufl.TestFunction(V_t), ufl.TrialFunction(V_t)
    interpolate_test = interpolated in ("test", "both")
    interpolate_trial = interpolated in ("trial", "both")

    def form(u, v):
        """A mass and a stiffness term, reading values and derivatives."""
        return (ufl.inner(u, v) + ufl.inner(ufl.grad(u), ufl.grad(v))) * ufl.dx

    a = form(
        ufl.Interpolate(u, V_t) if interpolate_trial else u,
        ufl.Interpolate(v, V_t) if interpolate_test else v,
    )
    a_ref = form(u_t if interpolate_trial else u, v_t if interpolate_test else v)

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [a, a_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    cell_type = module.lib.cell

    def tabulate(form, shape):
        offsets = form.form_integral_offsets
        assert offsets[cell_type + 1] - offsets[cell_type] == 1
        integral = form.form_integrals[offsets[cell_type]]
        A = np.zeros(shape, dtype=dtype)
        getattr(integral, f"tabulate_tensor_{dtype}")(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    def interpolation_operator(element):
        """The operator into the target space, on the dofs of the blocked element."""

        def base_element(element):
            return element.sub_elements[0] if element.block_size > 1 else element

        interpolation = basix.compute_interpolation_operator(
            base_element(element).basix_element, base_element(element_t).basix_element
        )
        return np.kron(interpolation, np.eye(element.block_size))

    M_0 = interpolation_operator(element_0) if interpolate_test else np.eye(element_0.dim)
    M_1 = interpolation_operator(element_1) if interpolate_trial else np.eye(element_1.dim)
    A = tabulate(compiled_forms[0], (element_0.dim, element_1.dim))
    A_ref = tabulate(compiled_forms[1], (M_0.shape[0], M_1.shape[0]))

    tol = np.finfo(dtype).eps * 100
    assert np.abs(A).max() > tol
    np.testing.assert_allclose(A, M_0.T @ A_ref @ M_1, rtol=tol, atol=tol)


@pytest.mark.parametrize("degrees", [(2, 1), (1, 2), (3, 2)])
@pytest.mark.parametrize("rank", [1, 2])
def test_interpolate_argument_tensor_product(compile_args, degrees, rank):
    """Interpolating an argument between tensor product spaces on a quadrilateral.

    Sum factorisation is on by default for these elements. The interpolated
    table includes the interpolation operator, which the 1D tables of the
    target element leave out, so it must not be sum factorised. The bilinear
    form tests against the target space, whose own table is sum factorised.
    """
    dtype = "float64"

    def element(degree, shape=None):
        e = basix.ufl.wrap_element(
            basix.create_tp_element(
                basix.ElementFamily.P,
                basix.CellType.quadrilateral,
                degree,
                basix.LagrangeVariant.gll_warped,
            )
        )
        return e if shape is None else basix.ufl.blocked_element(e, shape=shape)

    domain = ufl.Mesh(element(1, (2,)))
    element_a, element_t = element(degrees[0]), element(degrees[1])
    V_a = ufl.FunctionSpace(domain, element_a)
    V_t = ufl.FunctionSpace(domain, element_t)
    u_a, u_t, v_t = ufl.TrialFunction(V_a), ufl.TrialFunction(V_t), ufl.TestFunction(V_t)
    if rank == 1:
        J, J_ref = ufl.Interpolate(u_a, V_t) * ufl.dx, u_t * ufl.dx
    else:
        J, J_ref = v_t * ufl.Interpolate(u_a, V_t) * ufl.dx, v_t * u_t * ufl.dx

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    coords = np.array(
        [[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0], [2.3, 1.9, 0.0]], dtype=np.float64
    ).flatten()
    ffi = module.ffi

    def tabulate(form, shape):
        offsets = form.form_integral_offsets
        integral = form.form_integrals[offsets[module.lib.cell]]
        A = np.zeros(shape, dtype=dtype)
        getattr(integral, f"tabulate_tensor_{dtype}")(
            ffi.cast(f"{dtype_to_c_type(dtype)} *", A.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.cast("double *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    # The axis of the test function, for the bilinear form, then of the trial function
    test_axis = (element_t.dim,) if rank == 2 else ()
    A = tabulate(compiled_forms[0], (*test_axis, element_a.dim))
    A_ref = tabulate(compiled_forms[1], (*test_axis, element_t.dim))
    interpolation = basix.compute_interpolation_operator(
        element_a.basix_element, element_t.basix_element
    )
    np.testing.assert_allclose(A, A_ref @ interpolation, atol=1e-13)


def test_interpolate_derivative_with_test_function(compile_args):
    """Differentiate through an interpolation in a form that has a test function.

    The derivative pushes the trial function into the interpolation, whose
    Coargument has the same number as the form's test function. Interpolating
    ``w`` into its own space is the identity, so the Jacobian must be that of
    ``w * w * v * dx``.
    """
    dtype = "float64"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    w, v, du = ufl.Coefficient(V), ufl.TestFunction(V), ufl.TrialFunction(V)
    J = ufl.derivative(w * ufl.Interpolate(w, V) * v * ufl.dx, w, du)
    J_ref = ufl.derivative(w * w * v * ufl.dx, w, du)

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.3, 0.0], [0.1, 1.5, 0.0]]).flatten()
    w_values = np.random.default_rng(0).random(V.ufl_element().dim)
    constants = np.zeros(1)
    ffi = module.ffi

    def tabulate(form):
        offsets = form.form_integral_offsets
        integral = form.form_integrals[offsets[module.lib.cell]]
        A = np.zeros((V.ufl_element().dim,) * 2, dtype=dtype)
        integral.tabulate_tensor_float64(
            ffi.cast("double *", A.ctypes.data),
            ffi.cast("double *", w_values.ctypes.data),
            ffi.cast("double *", constants.ctypes.data),
            ffi.cast("double *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    A, A_ref = tabulate(compiled_forms[0]), tabulate(compiled_forms[1])
    assert np.abs(A_ref).max() > 1e-8
    np.testing.assert_allclose(A, A_ref, atol=1e-13)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
def test_interpolate_derivative(compile_args, dtype):
    """Differentiating an interpolation is the interpolation of the derivative.

    Interpolating into the space a Lagrange coefficient already lives in is the
    identity, so the differentiated interpolation must reproduce the plain
    derivative exactly.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element = basix.ufl.element("Lagrange", cell, 1, shape=(2,))
    V = ufl.FunctionSpace(domain, element)

    u = ufl.Coefficient(V)
    x = ufl.SpatialCoordinate(domain)
    expression = x + u

    J = ufl.derivative(ufl.Interpolate(expression, V), u)[0] * ufl.dx
    J_ref = ufl.derivative(expression, u)[0] * ufl.dx

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    w = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7], dtype=dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    cell_type = module.lib.cell

    tensors = []
    for form in compiled_forms:
        offsets = form.form_integral_offsets
        assert offsets[cell_type + 1] - offsets[cell_type] == 1
        integral = form.form_integrals[offsets[cell_type]]
        A = np.zeros(element.dim, dtype=dtype)
        kernel = getattr(integral, f"tabulate_tensor_{dtype}")
        kernel(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.cast(f"{c_type} *", w.ctypes.data),
            ffi.NULL,
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        tensors.append(A)

    tol = np.finfo(dtype).eps * 100
    assert np.abs(tensors[1]).max() > tol
    np.testing.assert_allclose(tensors[0], tensors[1], atol=tol)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("ngrads", [0, 1])
@pytest.mark.parametrize("rank", [1, 2])
def test_interpolate_derivative_in_form(compile_args, dtype, ngrads, rank):
    """Differentiating a form holding an interpolation applies the chain rule.

    ``d/du F(I(u))[v] = <dF/dI, I(v)>``, and interpolation is linear, so the
    tensor must be the reference interpolation operator applied to the tensor of
    the same form posed on the target space. A functional is differentiated
    along a test function, and a residual along a trial function, giving its
    Jacobian.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element_a = basix.ufl.element("Lagrange", cell, 1)
    element_t = basix.ufl.element("Lagrange", cell, 2)
    V_a = ufl.FunctionSpace(domain, element_a)
    V_t = ufl.FunctionSpace(domain, element_t)

    u, v = ufl.Coefficient(V_a), ufl.TestFunction(V_a)
    argument = ufl.TestFunction if rank == 1 else ufl.TrialFunction
    du, du_t = argument(V_a), argument(V_t)

    def differentiate(expression):
        """Take ``ngrads`` gradients and index back down to a scalar."""
        for _ in range(ngrads):
            expression = ufl.grad(expression)
        return expression[(0,) * ngrads] if ngrads else expression

    def form(expression):
        """The functional of the expression, or the residual testing it against ``v``."""
        return expression * ufl.dx if rank == 1 else expression * v * ufl.dx

    J = ufl.derivative(form(differentiate(ufl.Interpolate(u, V_t))), u, du)
    J_ref = form(differentiate(du_t))

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    cell_type = module.lib.cell

    def tabulate(form, shape):
        offsets = form.form_integral_offsets
        assert offsets[cell_type + 1] - offsets[cell_type] == 1
        integral = form.form_integrals[offsets[cell_type]]
        A = np.zeros(shape, dtype=dtype)
        getattr(integral, f"tabulate_tensor_{dtype}")(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    # The axis of the test function, for the Jacobian, then of the direction
    test_axis = (element_a.dim,) if rank == 2 else ()
    A = tabulate(compiled_forms[0], (*test_axis, element_a.dim))
    A_ref = tabulate(compiled_forms[1], (*test_axis, element_t.dim))

    interpolation = basix.compute_interpolation_operator(
        element_a.basix_element, element_t.basix_element
    )
    tol = np.finfo(dtype).eps * 100
    assert np.abs(A).max() > tol
    np.testing.assert_allclose(A, A_ref @ interpolation, atol=tol)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("derivative", [False, True], ids=["value", "curl-div"])
@pytest.mark.parametrize("rank", [1, 2])
@pytest.mark.parametrize("family", ["N1curl", "RT"])
def test_interpolate_derivative_piola(compile_args, dtype, derivative, rank, family):
    """Differentiating a form holding an interpolation into a Piola mapped space.

    As in `test_interpolate_derivative_in_form`, the derivative of ``F(I(u))``
    along ``du`` is ``F`` applied to ``I(du)``, the interpolation of the
    argument itself, so its table is known at compile time. The source and
    target are pulled back alike, so the Jacobian factors of their Piola maps
    cancel and the tensor is the reference interpolation operator applied to
    the tensor of the same form posed on the target space. The form reads a
    component of the interpolation, or the natural derivative of the space: the
    curl of an N1curl function and the divergence of an RT function.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element_a = basix.ufl.element(family, cell, 1)
    element_t = basix.ufl.element(family, cell, 2)
    V_a = ufl.FunctionSpace(domain, element_a)
    V_t = ufl.FunctionSpace(domain, element_t)

    u, v = ufl.Coefficient(V_a), ufl.TestFunction(V_a)
    argument = ufl.TestFunction if rank == 1 else ufl.TrialFunction
    du, du_t = argument(V_a), argument(V_t)
    operator = ufl.curl if family == "N1curl" else ufl.div

    def scalarise(f):
        """The natural derivative of the space, or a component of the value."""
        return operator(f) if derivative else f[0]

    def form(f):
        """The functional of ``f``, or the residual testing it against ``v``."""
        return f * ufl.dx if rank == 1 else f * scalarise(v) * ufl.dx

    J = ufl.derivative(form(scalarise(ufl.Interpolate(u, V_t))), u, du)
    J_ref = form(scalarise(du_t))

    # The interpolation must be recognised as that of the argument itself
    assert not ffcx.analysis.analyze_ufl_objects([J], np.dtype(dtype)).sub_kernels

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    cell_type = module.lib.cell

    def tabulate(form, shape):
        offsets = form.form_integral_offsets
        assert offsets[cell_type + 1] - offsets[cell_type] == 1
        integral = form.form_integrals[offsets[cell_type]]
        A = np.zeros(shape, dtype=dtype)
        getattr(integral, f"tabulate_tensor_{dtype}")(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    # The axis of the test function, for the Jacobian, then of the direction
    test_axis = (element_a.dim,) if rank == 2 else ()
    A = tabulate(compiled_forms[0], (*test_axis, element_a.dim))
    A_ref = tabulate(compiled_forms[1], (*test_axis, element_t.dim))

    interpolation = basix.compute_interpolation_operator(
        element_a.basix_element, element_t.basix_element
    )
    tol = np.finfo(dtype).eps * 100
    assert np.abs(A).max() > tol
    np.testing.assert_allclose(A, A_ref @ interpolation, atol=tol)


@pytest.mark.parametrize("family", ["N1curl", "RT"])
def test_interpolate_derivative_piola_runtime_table(family):
    """A table built for each cell is not mapped onto a Piola mapped target.

    Differentiating a non-linear interpolated expression gives an interpolation
    that is only linear in the argument, whose table is built for each cell.
    That is only supported for an identity mapped target.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    V_a = ufl.FunctionSpace(domain, basix.ufl.element(family, cell, 1))
    V_t = ufl.FunctionSpace(domain, basix.ufl.element(family, cell, 2))
    u, du = ufl.Coefficient(V_a), ufl.TestFunction(V_a)
    J = ufl.derivative(ufl.Interpolate(ufl.inner(u, u) * u, V_t)[0] * ufl.dx, u, du)
    with pytest.raises(NotImplementedError, match="only supported when the expression is"):
        ffcx.analysis.analyze_ufl_objects([J], np.float64)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("ngrads", [0, 1])
@pytest.mark.parametrize("block_size", [1, 2])
def test_interpolate_derivative_runtime_table(compile_args, dtype, ngrads, block_size):
    """A non-linear interpolated expression needs a table built for each cell.

    Differentiating it gives an interpolation of an expression that is linear in
    the argument but scaled by coefficients, so the argument's table depends on
    the cell. Check it against the interpolation applied by hand:
    ``A[j] = sum_I l_I(dh/du[phi_j]) * ref[I]``, with ``l_I`` the target
    element's interpolation functionals and ``ref`` the same form on the target.
    """
    cell = "triangle"
    shape = () if block_size == 1 else (block_size,)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element_a = basix.ufl.element("Lagrange", cell, 1, shape=shape)
    # A higher degree target only makes sense when it is not blocked, since a
    # blocked target must keep the value shape of the interpolated expression.
    element_t = basix.ufl.element("Lagrange", cell, 2 if block_size == 1 else 1, shape=shape)
    V_a = ufl.FunctionSpace(domain, element_a)
    V_t = ufl.FunctionSpace(domain, element_t)

    u, du = ufl.Coefficient(V_a), ufl.TestFunction(V_a)
    v_t = ufl.TrialFunction(V_t)

    # A non-linear expression, so that d/du brings `u` into the interpolation.
    if block_size == 1:
        expression = u * u
    else:
        expression = ufl.as_vector([u[i] * u[i] for i in range(block_size)])

    def scalarise(f):
        """Take ``ngrads`` gradients and index down to a scalar."""
        if block_size > 1:
            f = f[0]
        for _ in range(ngrads):
            f = ufl.grad(f)
        return f[(0,) * ngrads] if ngrads else f

    J = ufl.derivative(scalarise(ufl.Interpolate(expression, V_t)) * ufl.dx, u, du)
    J_ref = scalarise(v_t) * ufl.dx

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    w = np.arange(1, element_a.dim + 1, dtype=dtype) / 2
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    cell_type = module.lib.cell

    def tabulate(form, size, coefficients=None):
        integral = form.form_integrals[form.form_integral_offsets[cell_type]]
        A = np.zeros(size, dtype=dtype)
        getattr(integral, f"tabulate_tensor_{dtype}")(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.NULL if coefficients is None else ffi.cast(f"{c_type} *", coefficients.ctypes.data),
            ffi.NULL,
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    A = tabulate(compiled_forms[0], element_a.dim, w)
    A_ref = tabulate(compiled_forms[1], element_t.dim)

    # Apply the interpolation by hand. Blocked elements interpolate each block
    # component with the scalar sub-element's functionals.
    def base_element(element):
        return element.sub_elements[0] if element.block_size > 1 else element

    scalar_t = base_element(element_t).basix_element
    scalar_a = base_element(element_a).basix_element
    matrix, points = scalar_t.interpolation_matrix, scalar_t.points
    phi = scalar_a.tabulate(0, points)[0][:, :, 0]
    u_at_points = np.stack([phi @ w[c::block_size] for c in range(block_size)], axis=1)

    # d(expression[c])/d(dof j), with j the blocked dof (node, component).
    derivative = np.zeros((points.shape[0], element_a.dim, block_size))
    for j in range(element_a.dim):
        node, component = divmod(j, block_size)
        derivative[:, j, component] = 2 * u_at_points[:, component] * phi[:, node]

    interpolated = np.zeros((element_t.dim, element_a.dim))
    for i in range(matrix.shape[0]):
        for c in range(block_size):
            interpolated[i * block_size + c] = matrix[i] @ derivative[:, :, c]

    tol = np.finfo(dtype).eps * 100
    assert np.abs(A).max() > tol
    np.testing.assert_allclose(A, interpolated.T @ A_ref, atol=tol, rtol=100 * tol)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("flat_component", [0, 1, 2])
def test_interpolate_derivative_mixed_shapes(compile_args, dtype, flat_component):
    """A runtime table may interpolate into a space of a different value shape.

    The proxy argument then stands on the target space, for its value shape, but
    its table is indexed by the degrees of freedom of the argument it replaces.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element_a = basix.ufl.element("Lagrange", cell, 1, shape=(2,))
    element_t = basix.ufl.element("Lagrange", cell, 1, shape=(3,))
    V_a = ufl.FunctionSpace(domain, element_a)
    V_t = ufl.FunctionSpace(domain, element_t)

    r, dr = ufl.Coefficient(V_a), ufl.TestFunction(V_a)
    v_t = ufl.TrialFunction(V_t)
    expression = ufl.as_vector([r[0] * r[0], r[1] * r[1], r[0] * r[1]])

    J = ufl.derivative(ufl.Interpolate(expression, V_t)[flat_component] * ufl.dx, r, dr)
    J_ref = v_t[flat_component] * ufl.dx

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    w = np.array([1.0, 2.0, 0.5, 1.5, 2.5, 0.25], dtype=dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    cell_type = module.lib.cell

    def tabulate(form, size, coefficients=None):
        integral = form.form_integrals[form.form_integral_offsets[cell_type]]
        A = np.zeros(size, dtype=dtype)
        getattr(integral, f"tabulate_tensor_{dtype}")(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.NULL if coefficients is None else ffi.cast(f"{c_type} *", coefficients.ctypes.data),
            ffi.NULL,
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    A = tabulate(compiled_forms[0], element_a.dim, w)
    A_ref = tabulate(compiled_forms[1], element_t.dim)

    # Apply the interpolation by hand.
    scalar_t = element_t.sub_elements[0].basix_element
    scalar_a = element_a.sub_elements[0].basix_element
    matrix, points = scalar_t.interpolation_matrix, scalar_t.points
    phi = scalar_a.tabulate(0, points)[0][:, :, 0]
    r_at_points = np.stack([phi @ w[c::2] for c in range(2)], axis=1)

    # d(expression[component])/d(dof j), with j the blocked dof (node, component).
    derivative = np.zeros((points.shape[0], element_a.dim, 3))
    for j in range(element_a.dim):
        node, component = divmod(j, 2)
        if component == 0:
            derivative[:, j, 0] = 2 * r_at_points[:, 0] * phi[:, node]
            derivative[:, j, 2] = r_at_points[:, 1] * phi[:, node]
        else:
            derivative[:, j, 1] = 2 * r_at_points[:, 1] * phi[:, node]
            derivative[:, j, 2] = r_at_points[:, 0] * phi[:, node]

    interpolated = np.zeros((element_t.dim, element_a.dim))
    for i in range(matrix.shape[0]):
        for c in range(3):
            interpolated[i * 3 + c] = matrix[i] @ derivative[:, :, c]

    tol = np.finfo(dtype).eps * 100
    assert np.abs(A).max() > tol
    np.testing.assert_allclose(A, interpolated.T @ A_ref, atol=tol, rtol=100 * tol)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
def test_expression_two_arguments(compile_args, dtype):
    """An expression may hold two arguments.

    The interpolation of an expression bilinear in two arguments needs one, and
    the tensor written is ``[point][component][dof][dof]``.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element = basix.ufl.element("Lagrange", cell, 1)
    points = basix.ufl.element("Lagrange", cell, 2).basix_element.points
    V = ufl.FunctionSpace(domain, element)
    v, w = ufl.TestFunction(V), ufl.TrialFunction(V)

    compiled, module, _code = ffcx.codegeneration.jit.compile_expressions(
        [(2 * v * w, points)],
        options={"scalar_type": dtype},
        cffi_extra_compile_args=compile_args,
    )
    expression = compiled[0]
    assert expression.num_points == points.shape[0]

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype).flatten()
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    A = np.zeros(points.shape[0] * element.dim**2, dtype=dtype)
    getattr(expression, f"tabulate_tensor_{dtype}")(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    phi = element.basix_element.tabulate(0, points)[0][:, :, 0]
    tol = np.finfo(dtype).eps * 100
    np.testing.assert_allclose(
        A.reshape(points.shape[0], element.dim, element.dim),
        2 * phi[:, :, None] * phi[:, None, :],
        atol=tol,
    )


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("ngrads", [0, 1])
@pytest.mark.parametrize("block_size", [1, 2])
def test_interpolate_second_derivative(compile_args, dtype, ngrads, block_size):
    """The Jacobian of an interpolated non-linear expression is a dense block.

    Differentiating twice gives an interpolation of an expression bilinear in
    both arguments, so it takes one element tensor axis per argument rather than
    being a product of a factor per axis. Check it against a finite difference
    of the residual, which is verified separately.
    """
    cell = "triangle"
    shape = () if block_size == 1 else (block_size,)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element_a = basix.ufl.element("Lagrange", cell, 1, shape=shape)
    element_t = basix.ufl.element("Lagrange", cell, 2 if block_size == 1 else 1, shape=shape)
    V_a = ufl.FunctionSpace(domain, element_a)
    V_t = ufl.FunctionSpace(domain, element_t)

    u = ufl.Coefficient(V_a)
    du, ddu = ufl.TestFunction(V_a), ufl.TrialFunction(V_a)

    if block_size == 1:
        expression = u * u
    else:
        expression = ufl.as_vector([u[i] * u[i] for i in range(block_size)])

    integrand = ufl.Interpolate(expression, V_t)
    if block_size > 1:
        integrand = integrand[0]
    for _ in range(ngrads):
        integrand = ufl.grad(integrand)
    if ngrads:
        integrand = integrand[(0,) * ngrads]

    residual = ufl.derivative(integrand * ufl.dx, u, du)
    jacobian = ufl.derivative(residual, u, ddu)

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [residual, jacobian],
        options={"scalar_type": dtype},
        cffi_extra_compile_args=compile_args,
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    cell_type = module.lib.cell
    n = element_a.dim

    def tabulate(form, size, coefficients):
        integral = form.form_integrals[form.form_integral_offsets[cell_type]]
        A = np.zeros(size, dtype=dtype)
        getattr(integral, f"tabulate_tensor_{dtype}")(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.cast(f"{c_type} *", coefficients.ctypes.data),
            ffi.NULL,
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    w = (np.arange(1, n + 1) / 3).astype(dtype)
    A = tabulate(compiled_forms[1], n * n, w).reshape(n, n)

    # Finite difference of the residual, column by column.
    step = 1e-6 if dtype == "float64" else 1e-3
    expected = np.zeros((n, n))
    for k in range(n):
        direction = np.zeros(n, dtype=dtype)
        direction[k] = step
        expected[:, k] = (
            tabulate(compiled_forms[0], n, w + direction)
            - tabulate(compiled_forms[0], n, w - direction)
        ) / (2 * step)

    scale = max(np.abs(expected).max(), 1e-30)
    assert np.abs(A).max() > 1e-10
    assert np.abs(A - expected).max() / scale < (1e-6 if dtype == "float64" else 1e-2)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("degrees", [(1, 2), (2, 1)], ids=["P1xP2", "P2xP1"])
def test_interpolate_second_derivative_mixed_spaces(compile_args, dtype, degrees):
    """The dense block of a Hessian takes each argument's own dof extent.

    Differentiating ``Interpolate(u*u, P2)*dx`` twice, with the test and trial
    functions in different spaces, gives one table spanning both element
    tensor axes. Its extent on each axis is that of the argument the axis
    belongs to, so the two have to be told apart when the block is written
    into the element tensor. The Hessian is the interpolation of ``2 du ddu``
    into P2, so for Lagrange it is the values of ``2 du ddu`` at the P2 nodes
    weighted by the integrals of the P2 basis functions.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    spaces = {
        degree: ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", cell, degree))
        for degree in (1, 2)
    }
    V_test, V_trial = (spaces[degree] for degree in degrees)
    P1, P2 = spaces[1], spaces[2]

    u = ufl.Coefficient(P1)
    du, ddu = ufl.TestFunction(V_test), ufl.TrialFunction(V_trial)
    residual = ufl.derivative(ufl.Interpolate(u * u, P2) * ufl.dx, u, du)
    jacobian = ufl.derivative(residual, u, ddu)

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [jacobian], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    form = compiled_forms[0]
    integral = form.form_integrals[form.form_integral_offsets[module.lib.cell]]
    shape = (V_test.ufl_element().dim, V_trial.ufl_element().dim)
    size = shape[0] * shape[1]
    w = (np.arange(1, 4) / 3).astype(dtype)
    # Padded, so that a kernel writing past the element tensor is caught
    A = np.zeros(2 * size, dtype=dtype)
    getattr(integral, f"tabulate_tensor_{dtype}")(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.NULL,
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    assert not np.any(A[size:])

    target = P2.ufl_element().basix_element
    phi_test = V_test.ufl_element().basix_element.tabulate(0, target.points)[0, :, :, 0]
    phi_trial = V_trial.ufl_element().basix_element.tabulate(0, target.points)[0, :, :, 0]
    points, weights = basix.make_quadrature(basix.CellType.triangle, 2)
    detJ = abs(np.linalg.det(np.column_stack([coords[1] - coords[0], coords[2] - coords[0]])[:2]))
    integrals = detJ * (weights @ target.tabulate(0, points)[0, :, :, 0])
    expected = 2 * np.einsum(
        "p,pq,qi,qj->ij", integrals, target.interpolation_matrix, phi_test, phi_trial
    )
    tol = np.finfo(dtype).eps * 100
    np.testing.assert_allclose(A[:size].reshape(shape), expected, atol=tol)


@pytest.mark.parametrize(
    "target", [("Lagrange", 1), ("Lagrange", 2), ("N1curl", 1)], ids=["P1", "P2", "N1curl1"]
)
def test_shape_derivative_of_interpolation(compile_args, target):
    """The shape derivative of a functional of an interpolation.

    A simple shape optimisation problem: the misfit ``J = 1/2 |I_W u - u_d|^2 dx``
    of a given state ``u`` interpolated into ``W`` and a target ``u_d`` given in
    terms of ``x``. Its shape derivative ``dJ[V] = derivative(J, x, V)`` is
    checked by a Taylor test on one cell: moving the vertices by ``t * V`` with
    the dofs of ``u`` kept, ``J`` changes by ``t * dJ[V]`` up to ``O(t**2)``.
    The dofs of the interpolation move with the cell, through its interpolation
    points and, for a Piola mapped target, its Jacobian. For the P2 target the
    interpolation is the identity, so ``dJ`` must also equal the derivative of
    the functional without it.
    """
    family, degree = target
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    element = basix.ufl.element("Lagrange", "triangle", 2)
    u = ufl.Coefficient(ufl.FunctionSpace(domain, element))
    W = ufl.FunctionSpace(domain, basix.ufl.element(family, "triangle", degree))
    x = ufl.SpatialCoordinate(domain)
    # The shape derivative is taken along the test function of the coordinate space
    V = ufl.TestFunction(ufl.FunctionSpace(domain, domain.ufl_coordinate_element()))
    if family == "N1curl":
        state, u_d = ufl.grad(u), ufl.as_vector([x[1], -x[0] * x[1]])
    else:
        state, u_d = u, x[0] * x[1]

    def functional(q):
        return 0.5 * ufl.inner(q - u_d, q - u_d) * ufl.dx

    J = functional(ufl.Interpolate(state, W))
    forms = [J, ufl.derivative(J, x, V), ufl.derivative(functional(state), x, V)]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    ffi = module.ffi
    rng = np.random.default_rng(3)
    u_dofs = rng.uniform(0.5, 1.5, element.dim)
    vertices = np.array([[0.1, 0.2], [2.0, 0.3], [0.4, 1.7]])
    direction = rng.uniform(-1, 1, vertices.shape)

    def tabulate(form, t, size):
        """Tabulate on the cell with its vertices moved by ``t`` along the direction."""
        coords = np.pad(vertices + t * direction, ((0, 0), (0, 1))).flatten()
        A = np.zeros(size)
        _integral(module, form).tabulate_tensor_float64(
            ffi.cast("double *", A.ctypes.data),
            ffi.cast("double *", u_dofs.ctypes.data),
            ffi.NULL,
            ffi.cast("double *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    # The dofs of the coordinate space are the vertex displacements, interleaved
    num_dofs = direction.size
    dJ = tabulate(compiled_forms[1], 0.0, num_dofs) @ direction.flatten()
    assert abs(dJ) > 1e-3
    J0 = tabulate(compiled_forms[0], 0.0, 1)[0]
    steps = 1e-2 / 2 ** np.arange(4)
    residuals = np.array([abs(tabulate(compiled_forms[0], t, 1)[0] - J0 - t * dJ) for t in steps])
    assert residuals[0] > 1e-8
    rates = np.log2(residuals[:-1] / residuals[1:])
    np.testing.assert_allclose(rates, 2, atol=0.05)

    if target == ("Lagrange", 2):
        dJ_dofs, dJ_ref = (tabulate(form, 0.0, num_dofs) for form in compiled_forms[1:])
        np.testing.assert_allclose(dJ_dofs, dJ_ref, atol=1e-13)


@pytest.mark.parametrize("do_cancel_jacobian_products", [False, True])
@pytest.mark.parametrize(
    "scalar_type, complex_node",
    [(np.float64, None), (np.float64, ufl.conj), (np.float64, ufl.real), (np.complex128, None)],
    ids=["float64", "float64-conj", "float64-real", "complex128"],
)
@pytest.mark.parametrize("measure", ["dx", "ds", "expression"])
@pytest.mark.parametrize("family", ["Lagrange", "N1curl", "RT", "Regge", "HHJ"])
def test_interpolate_argument_is_compile_time(
    family, measure, scalar_type, complex_node, do_cancel_jacobian_products
):
    """Interpolating an argument itself must not need a table built per cell.

    UFL leaves the argument mapped by its own pull back and the inverse of the
    target's, with the geometry of a Piola map lowered, and reduces it to its
    reference value only if Jacobian products are cancelled. The analysis
    recognises it either way, and rewrites it into its reference value.

    For a real scalar type UFL also removes a ``conj`` or ``real`` around the
    argument, so the operand only becomes the argument during processing. It is
    recognised all the same, which a Piola mapped target relies on: only the
    argument itself can be interpolated into one.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    V_a = ufl.FunctionSpace(domain, basix.ufl.element(family, cell, 2))
    V_t = ufl.FunctionSpace(domain, basix.ufl.element(family, cell, 1))
    argument = ufl.TestFunction(V_a)

    operand = argument if complex_node is None else complex_node(argument)
    interpolated = ufl.Interpolate(operand, V_t)
    integrand = sum(interpolated[i] for i in np.ndindex(interpolated.ufl_shape))
    if measure == "expression":
        points = np.array([[0.25, 0.25]])
        analysis = ffcx.analysis.analyze_ufl_objects([(integrand, points)], scalar_type)
        integrands = [expression for expression, _, _ in analysis.expressions]
    else:
        form = ufl.conj(integrand) * ufl.Measure(measure, domain=domain)
        analysis = ffcx.analysis.analyze_ufl_objects(
            [form], scalar_type, do_cancel_jacobian_products
        )
        integrands = [
            integral.integrand()
            for form_data in analysis.form_data
            for integral_data in form_data.integral_data
            for integral in integral_data.integrals
        ]

    interpolations = {
        node
        for integrand in integrands
        for node in ufl.corealg.traversal.unique_pre_traversal(integrand)
        if isinstance(node, ufl.Interpolate)
    }
    assert len(interpolations) == 1
    (node,) = interpolations
    assert node.ufl_operands[0] == ufl.classes.ReferenceValue(argument)
    assert ffcx.analysis.interpolated_argument(node) == argument
    assert not ffcx.analysis.interpolation_has_runtime_table(node)
    assert not analysis.sub_kernels


def _integral(module, form, subdomain_id=-1):
    """Get the cell integral of a compiled form over a subdomain."""
    cell = module.lib.cell
    offsets = form.form_integral_offsets
    for i in range(offsets[cell], offsets[cell + 1]):
        if form.form_integral_ids[i] == subdomain_id:
            return form.form_integrals[i]
    raise KeyError(subdomain_id)


def _tabulate_scalar(module, integral, coefficients=None, constants=None, dtype="float64"):
    """Tabulate a functional on an affine triangle."""
    ffi = module.ffi
    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]], dtype=xdtype).flatten()
    A = np.zeros(1, dtype=dtype)

    def pointer(array, c_type):
        return ffi.NULL if array is None else ffi.cast(f"{c_type} *", array.ctypes.data)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    getattr(integral, f"tabulate_tensor_{dtype}")(
        pointer(A, c_type),
        pointer(coefficients, c_type),
        pointer(constants, c_type),
        pointer(coords, c_xtype),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    return A[0]


def test_interpolation_shared_between_forms(compile_args):
    """An interpolation in two forms, or two that differ only in coefficients.

    The kernel of the interpolation is generated once, and each interpolation
    is evaluated with its own coefficients. Interpolating into the space the
    coefficients live in is the identity.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    u, w = ufl.Coefficient(V), ufl.Coefficient(V)
    forms = [
        (ufl.Interpolate(u, V) + 2 * ufl.Interpolate(w, V)) * ufl.dx,
        3 * ufl.Interpolate(u, V) * ufl.dx,
        (u + 2 * w) * ufl.dx,
    ]
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    rng = np.random.default_rng(1)
    u_values, w_values = rng.random(6), rng.random(6)
    both = np.concatenate([u_values, w_values])
    A = _tabulate_scalar(module, _integral(module, compiled[0]), both)
    B = _tabulate_scalar(module, _integral(module, compiled[1]), u_values)
    A_ref = _tabulate_scalar(module, _integral(module, compiled[2]), both)
    np.testing.assert_allclose(A, A_ref, rtol=1e-13)
    w_ref = _tabulate_scalar(module, _integral(module, compiled[2]), np.r_[0 * u_values, w_values])
    np.testing.assert_allclose(B, 3 * (A_ref - w_ref), rtol=1e-13)


def test_interpolation_constants_in_form(compile_args):
    """The constants of an interpolated expression are passed in its own order.

    The form holds a constant that the interpolation does not, ahead of those
    that it does.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    u = ufl.Coefficient(V)
    c0, c1, c2 = ufl.Constant(domain), ufl.Constant(domain, shape=(2,)), ufl.Constant(domain)
    forms = [
        c0 * ufl.Interpolate(c2 * u + c1[1], V) * ufl.dx,
        c0 * (c2 * u + c1[1]) * ufl.dx,
    ]
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    u_values = np.random.default_rng(2).random(6)
    constants = np.array([2.0, 3.0, -0.5, 0.25])
    A, A_ref = (
        _tabulate_scalar(module, _integral(module, form), u_values, constants) for form in compiled
    )
    np.testing.assert_allclose(A, A_ref, rtol=1e-13)


def test_interpolation_in_one_of_several_integrals(compile_args):
    """Interpolations in different integrals of one form.

    Each integral evaluates the proxy coefficients it holds, where it reads
    them. The interpolation in the last integral reads a coefficient that the
    form holds elsewhere, which must be passed to that integral too.
    Interpolating into the space a coefficient lives in is the identity.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V1 = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 1))
    V2 = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    u, z = ufl.Coefficient(V2), ufl.Coefficient(V1)
    forms = [
        z * ufl.dx(1) + ufl.Interpolate(u, V2) * ufl.dx(2) + ufl.Interpolate(z, V1) * ufl.dx(3),
        (u + z) * ufl.dx,
    ]
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    form = compiled[0]
    assert form.num_coefficients == 2
    rng = np.random.default_rng(3)
    u_values, z_values = rng.random(6), rng.random(3)
    coefficients = np.concatenate([u_values, z_values])

    def reference(u, z):
        return _tabulate_scalar(module, _integral(module, compiled[1]), np.concatenate([u, z]))

    for subdomain, expected, enabled in [
        (1, reference(0 * u_values, z_values), [False, True]),
        (2, reference(u_values, 0 * z_values), [True, False]),
        (3, reference(0 * u_values, z_values), [False, True]),
    ]:
        integral = _integral(module, form, subdomain)
        assert [integral.enabled_coefficients[i] for i in range(2)] == enabled
        A = _tabulate_scalar(module, integral, coefficients)
        np.testing.assert_allclose(A, expected, rtol=1e-13)


def test_two_interpolations_in_one_integral(compile_args):
    """Each proxy coefficient of an integral is evaluated, not only the first."""
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V1 = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 1))
    V2 = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    u, z = ufl.Coefficient(V2), ufl.Coefficient(V1)
    forms = [
        (ufl.Interpolate(u, V2) - 2 * ufl.Interpolate(z, V1)) * ufl.dx,
        (u - 2 * z) * ufl.dx,
    ]
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    rng = np.random.default_rng(5)
    coefficients = np.concatenate([rng.random(6), rng.random(3)])
    A, A_ref = (
        _tabulate_scalar(module, _integral(module, form), coefficients) for form in compiled
    )
    np.testing.assert_allclose(A, A_ref, rtol=1e-13)


def test_nested_interpolation_in_form(compile_args):
    """An interpolation of an expression holding another interpolation.

    Both interpolate into the space the coefficient lives in, so they are the
    identity.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V1 = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 1))
    V2 = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    z = ufl.Coefficient(V1)
    forms = [
        ufl.Interpolate(ufl.Interpolate(z, V1) * z, V2) * ufl.dx,
        z * z * ufl.dx,
    ]
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    z_values = np.random.default_rng(4).random(3)
    A, A_ref = (_tabulate_scalar(module, _integral(module, form), z_values) for form in compiled)
    np.testing.assert_allclose(A, A_ref, rtol=1e-13)


def test_interpolation_on_interior_facet(compile_args):
    """Interpolations restricted to either side of an interior facet.

    Interpolating into the space the argument or coefficient lives in is the
    identity, so each form must match the one without the interpolation. A
    bilinear form restricts its two interpolations to each combination of
    sides. The dofs of a proxy coefficient are built on each of the two cells.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    v, u, f = ufl.TestFunction(V), ufl.TrialFunction(V), ufl.Coefficient(V)
    x = ufl.SpatialCoordinate(domain)

    # Two triangles sharing the facet (1, 0) -> (0, 1), local facet 0 of both,
    # with the same vertex order so neither side is permuted.
    coords = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64
    ).flatten()
    entity = np.array([0, 0], dtype=np.int32)
    perm = np.array([0, 0], dtype=np.uint8)
    # Far apart on the two sides, so reading the wrong side shows
    w = np.concatenate([np.random.default_rng(5).random(6), 10 + np.arange(6.0)])

    def tabulate(form, shape):
        (compiled,), module, _ = ffcx.codegeneration.jit.compile_forms(
            [form], cffi_extra_compile_args=compile_args
        )
        ffi = module.ffi
        offsets = compiled.form_integral_offsets
        integral = compiled.form_integrals[offsets[module.lib.interior_facet]]
        A = np.zeros(shape)
        integral.tabulate_tensor_float64(
            ffi.cast("double *", A.ctypes.data),
            ffi.cast("double *", w.ctypes.data),
            ffi.NULL,
            ffi.cast("double *", coords.ctypes.data),
            ffi.cast("int *", entity.ctypes.data),
            ffi.cast("uint8_t *", perm.ctypes.data),
            ffi.NULL,
        )
        return A

    for side in ("+", "-"):
        # The table of an interpolated argument, into that side's block
        A = tabulate(ufl.Interpolate(v, V)(side) * (1 + x[0]) * ufl.dS, 2 * 6)
        A_ref = tabulate(v(side) * (1 + x[0]) * ufl.dS, 2 * 6)
        assert np.abs(A_ref).max() > 1e-3
        np.testing.assert_allclose(A, A_ref, atol=1e-14)

    for side_v, side_u in itertools.product(("+", "-"), repeat=2):
        # The tables of both interpolated arguments, into their sides' blocks
        I_v, I_u = ufl.Interpolate(v, V)(side_v), ufl.Interpolate(u, V)(side_u)
        A = tabulate(I_v * I_u * (1 + x[0]) * ufl.dS, (2 * 6, 2 * 6))
        A_ref = tabulate(v(side_v) * u(side_u) * (1 + x[0]) * ufl.dS, (2 * 6, 2 * 6))
        assert np.abs(A_ref).max() > 1e-3
        np.testing.assert_allclose(A, A_ref, atol=1e-14)

    for side in ("+", "-"):
        # A proxy coefficient, with its dofs built on that side's cell
        A = tabulate(ufl.Interpolate(f, V)(side) * (1 + x[0]) * ufl.dS, 1)
        A_ref = tabulate(f(side) * (1 + x[0]) * ufl.dS, 1)
        np.testing.assert_allclose(A, A_ref, rtol=1e-13)

    # Both sides in one integrand, with the coefficients of the two cells apart
    A = tabulate(ufl.Interpolate(2 * f, V)("+") * ufl.Interpolate(f, V)("-") * ufl.dS, 1)
    A_ref = tabulate(2 * f("+") * f("-") * ufl.dS, 1)
    np.testing.assert_allclose(A, A_ref, rtol=1e-13)

    # A proxy coefficient in a space on a submesh of facets, built on each
    # side's cell. Its nodes are those of V on the facet, so it is the trace of
    # the interpolation into V. The expression reads the geometry, so building
    # the '-' side from the coordinates of the wrong cell shows.
    facets = ufl.Mesh(basix.ufl.element("Lagrange", "interval", 1, shape=(2,)))
    V_facet = ufl.FunctionSpace(facets, basix.ufl.element("Lagrange", "interval", 2))
    expression = f * (1 + x[1])
    for side in ("+", "-"):
        A = tabulate(ufl.Interpolate(expression, V_facet)(side) * (1 + x[0]) * ufl.dS, 1)
        A_ref = tabulate(ufl.Interpolate(expression, V)(side) * (1 + x[0]) * ufl.dS, 1)
        np.testing.assert_allclose(A, A_ref, rtol=1e-13)


@pytest.mark.parametrize("measure", ["ds", "dP", "dr"])
def test_interpolation_on_cell_entities(compile_args, measure):
    """Interpolations in exterior facet, vertex and ridge integrals.

    Each integral is tabulated on every entity of its kind of a tetrahedron.
    Interpolating into the space the argument or coefficient lives in is the
    identity, so each form must match the one without the interpolation. The
    dofs of a proxy coefficient are built on the cell, and the proxy is then
    read on the entity, also when it is interpolated into another space: into
    P1, a P2 coefficient keeps its vertex dofs.
    """
    cell = "tetrahedron"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(3,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", cell, 2))
    V1 = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", cell, 1))
    v, u, f, g = ufl.TestFunction(V), ufl.TrialFunction(V), ufl.Coefficient(V), ufl.Coefficient(V1)
    x = ufl.SpatialCoordinate(domain)
    dm = ufl.Measure(measure, domain=domain)
    integral_type = {"ds": "exterior_facet", "dP": "vertex", "dr": "ridge"}[measure]
    num_entities = {"ds": 4, "dP": 4, "dr": 6}[measure]

    coords = np.array(
        [[0.1, 0.2, 0.0], [2.0, 0.3, 0.1], [0.4, 1.7, 0.2], [0.3, 0.4, 1.5]], dtype=np.float64
    ).flatten()
    w = np.random.default_rng(5).random(V.ufl_element().dim)
    perm = np.zeros(1, dtype=np.uint8)

    def tabulate(form, shape, coefficients):
        """Tabulate the integral of a form on each entity of the cell."""
        (compiled,), module, _ = ffcx.codegeneration.jit.compile_forms(
            [form], cffi_extra_compile_args=compile_args
        )
        ffi = module.ffi
        offsets = compiled.form_integral_offsets
        integral = compiled.form_integrals[offsets[getattr(module.lib, integral_type)]]
        tensors = np.zeros((num_entities, *np.atleast_1d(shape)))
        for entity, A in enumerate(tensors):
            entity_index = np.array([entity], dtype=np.int32)
            integral.tabulate_tensor_float64(
                ffi.cast("double *", A.ctypes.data),
                ffi.cast("double *", coefficients.ctypes.data),
                ffi.NULL,
                ffi.cast("double *", coords.ctypes.data),
                ffi.cast("int *", entity_index.ctypes.data),
                ffi.cast("uint8_t *", perm.ctypes.data),
                ffi.NULL,
            )
        return tensors

    def check(form, form_ref, shape, coefficients=w, coefficients_ref=w, **tolerance):
        A, A_ref = tabulate(form, shape, coefficients), tabulate(form_ref, shape, coefficients_ref)
        assert np.abs(A_ref).max() > 1e-3
        # The entities differ, so the entity index is read
        assert not np.allclose(A_ref, A_ref[0])
        np.testing.assert_allclose(A, A_ref, **tolerance)

    dim = V.ufl_element().dim
    # The table of an interpolated argument
    check(ufl.Interpolate(v, V) * (1 + x[0]) * dm, v * (1 + x[0]) * dm, dim, atol=1e-14)
    # The tables of both interpolated arguments
    check(
        ufl.Interpolate(v, V) * ufl.Interpolate(u, V) * (1 + x[0]) * dm,
        v * u * (1 + x[0]) * dm,
        (dim, dim),
        atol=1e-14,
    )
    # A proxy coefficient, with its dofs built on the cell
    check(ufl.Interpolate(f, V) * (1 + x[0]) * dm, f * (1 + x[0]) * dm, 1, rtol=1e-13)
    # A proxy coefficient in another space: the vertex dofs of f come first
    check(
        ufl.Interpolate(f, V1) * (1 + x[0]) * dm,
        g * (1 + x[0]) * dm,
        1,
        coefficients_ref=w[: V1.ufl_element().dim],
        rtol=1e-13,
    )
    if measure == "dP":
        # A coefficient on a submesh of vertices is not supported yet in a
        # vertex integral, interpolated or not.
        return
    # A proxy coefficient in a space on a submesh of the entities, built on the
    # entity of the cell. Its nodes are those of V on the entity, so it is the
    # trace of the interpolation into V.
    entity_cell = {"ds": "triangle", "dr": "interval"}[measure]
    submesh = ufl.Mesh(basix.ufl.element("Lagrange", entity_cell, 1, shape=(3,)))
    V_entity = ufl.FunctionSpace(submesh, basix.ufl.element("Lagrange", entity_cell, 2))
    expression = f * (1 + x[1])
    check(
        ufl.Interpolate(expression, V_entity) * (1 + x[0]) * dm,
        ufl.Interpolate(expression, V) * (1 + x[0]) * dm,
        1,
        rtol=1e-13,
    )


@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        pytest.param(
            "complex128",
            marks=pytest.mark.skipif(
                sys.platform.startswith("win32"),
                reason="win32 platform does not support C99 _Complex numbers",
            ),
        ),
    ],
)
@pytest.mark.parametrize(
    "target", [("Lagrange", 1), ("DG", 0), ("N1curl", 1)], ids=["P1", "DG0", "N1curl1"]
)
def test_proxy_coefficient_restricted_to_minus(compile_args, dtype, target):
    """A proxy coefficient on the '-' side of an interior facet, and on both sides.

    The reference is a coefficient of the target space whose dofs are the
    interpolation, built by hand on each cell. The two cells differ in shape,
    so a Piola mapped target shows if the '-' dofs are built on the wrong
    cell, and their coefficients differ, so it shows if they are read from the
    wrong side.
    """
    family, degree = target
    vector = family == "N1curl"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    element_f = basix.ufl.element("Lagrange", "triangle", 2, shape=(2,) if vector else None)
    element_t = basix.ufl.element(family, "triangle", degree)
    V = ufl.FunctionSpace(domain, element_f)
    W = ufl.FunctionSpace(domain, element_t)
    f, g, v = ufl.Coefficient(V), ufl.Coefficient(W), ufl.TestFunction(V)
    expression = ufl.as_vector([f[0] * f[1], f[0] ** 2]) if vector else f**2
    interpolation = ufl.Interpolate(expression, W)

    def forms(q):
        return [
            ufl.inner(q("-"), v("-")) * ufl.dS,
            ufl.inner(q("+") - 2 * q("-"), v("+") + v("-")) * ufl.dS,
        ]

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms(interpolation) + forms(g),
        options={"scalar_type": dtype},
        cffi_extra_compile_args=compile_args,
    )

    # Two triangles of different shape sharing the facet (2, 0.3) -> (0.4, 1.7),
    # local facet 0 of both, with the same vertex order.
    cells = [
        np.array([[0.1, 0.2], [2.0, 0.3], [0.4, 1.7]]),
        np.array([[2.1, 1.9], [2.0, 0.3], [0.4, 1.7]]),
    ]
    rng = np.random.default_rng(7)
    f_dofs = [rng.uniform(0.5, 1.5, element_f.dim) for _ in cells]
    if np.issubdtype(dtype, np.complexfloating):
        f_dofs = [w + 1j * rng.uniform(-0.3, 0.3, w.size) for w in f_dofs]

    def interpolate_by_hand(cell, w):
        """The dofs of the interpolation on a cell, from the reference element."""
        points = element_t.basix_element.points
        phi = element_f.basix_element.tabulate(0, points)[0][:, :, 0]
        if not vector:
            return element_t.basix_element.interpolation_matrix @ (phi @ w) ** 2
        f0, f1 = phi @ w[0::2], phi @ w[1::2]
        values = np.stack([f0 * f1, f0**2], axis=1)
        # The covariant Piola pull back of an affine cell, J^T v
        J = (cell[1:] - cell[0]).T
        return element_t.basix_element.interpolation_matrix @ (values @ J).T.flatten()

    g_dofs = [interpolate_by_hand(cell, w) for cell, w in zip(cells, f_dofs)]
    coords = np.hstack([np.pad(c, ((0, 0), (0, 1))).flatten() for c in cells])
    entity = np.zeros(2, dtype=np.int32)
    perm = np.zeros(2, dtype=np.uint8)
    c_type = dtype_to_c_type(dtype)
    ffi = module.ffi

    def tabulate(form, w):
        offsets = form.form_integral_offsets
        integral = form.form_integrals[offsets[module.lib.interior_facet]]
        A = np.zeros(2 * element_f.dim, dtype=dtype)
        w = np.concatenate(w).astype(dtype)
        getattr(integral, f"tabulate_tensor_{dtype}")(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.cast(f"{c_type} *", w.ctypes.data),
            ffi.NULL,
            ffi.cast("double *", coords.ctypes.data),
            ffi.cast("int *", entity.ctypes.data),
            ffi.cast("uint8_t *", perm.ctypes.data),
            ffi.NULL,
        )
        return A

    for proxy_form, reference_form in zip(compiled_forms[:2], compiled_forms[2:]):
        A = tabulate(proxy_form, f_dofs)
        A_ref = tabulate(reference_form, g_dofs)
        assert np.abs(A_ref).max() > 1e-3
        np.testing.assert_allclose(A, A_ref, rtol=1e-12, atol=1e-12 * np.abs(A_ref).max())


def test_proxy_coefficients_are_not_form_coefficients():
    """A proxy coefficient is internal to the kernel, not a coefficient of the form.

    The coefficients it interpolates take its place: they are numbered, named,
    hashed and enabled as any other coefficient of the form.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    Q = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 1))
    u, f, g = ufl.Coefficient(V), ufl.Coefficient(Q), ufl.Coefficient(Q)
    v = ufl.TestFunction(V)
    # `f` is only read by the interpolation, `g` only by the second integral
    form = ufl.Interpolate(u**2 * f, Q) * v * ufl.dx + g * v * ufl.dx(1)

    options = ffcx.options.get_options()
    analysis = ffcx.analysis.analyze_ufl_objects([form], options["scalar_type"])
    ir = ffcx.ir.representation.compute_ir(analysis, {}, "test", options, None)

    (form_ir,) = ir.forms
    assert form_ir.num_coefficients == 3
    assert len(form_ir.coefficient_names) == 3
    assert form_ir.original_coefficient_positions == [0, 1, 2]
    assert len(form_ir.finite_element_hashes) == form_ir.rank + 3

    subdomains = dict(zip(form_ir.integral_names["cell"], form_ir.subdomain_ids["cell"]))
    enabled = {
        subdomains[integral.expression.name]: integral.enabled_coefficients
        for integral in ir.integrals
    }
    assert enabled == {-1: [True, True, False], 1: [False, False, True]}
    for integral in ir.integrals:
        assert not any(
            isinstance(c, ffcx.analysis.ProxyCoefficient)
            for c in integral.expression.coefficient_numbering
        )
    (sub_kernels,) = (integral.sub_kernels for integral in ir.integrals if integral.sub_kernels)
    assert len(sub_kernels) == 1


@pytest.mark.parametrize(
    "symmetric",
    [
        False,
        # Basix reports the reference value of a symmetric element with the
        # shape of the block, which UFL's inverse pull back cannot fill, see
        # https://github.com/FEniCS/basix/issues/1124
        pytest.param(
            True,
            marks=pytest.mark.xfail(
                strict=True,
                raises=NotImplementedError,
                reason="Basix symmetric elements report a block-shaped reference value, "
                "https://github.com/FEniCS/basix/issues/1124",
            ),
        ),
    ],
)
@pytest.mark.parametrize("dtype", ["float64", "float32"])
def test_interpolate_into_mixed_space(compile_args, dtype, symmetric):
    """Interpolating an expression into a mixed space.

    The space is vector DG1 x (RT x quadrature), with symmetric tensor P1 added
    to the inner mixed element if ``symmetric``. It holds a nested mixed element,
    a blocked sub-element whose interpolation matrix is not the identity, a Piola
    mapped sub-element, and a blocked quadrature sub-element, whose rule the integral
    takes. As the reference interpolates into each sub-space, it holds several
    proxy coefficients in one integral.

    It is checked twice: with a quadratic field, which the sub-spaces do not
    reproduce, against interpolations into each of them, the full tensor space
    for the symmetric one; and with a field whose integral the interpolation
    keeps, linear or the product x y at the quadrature points, against the exact
    integral of the field, which involves no interpolation.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    x = ufl.SpatialCoordinate(domain)
    DG1 = basix.ufl.element(
        "DG", "triangle", 1, shape=(2,), lagrange_variant=basix.LagrangeVariant.legendre
    )
    assert not DG1.basix_element.interpolation_is_identity
    RT = basix.ufl.element("RT", "triangle", 1)
    Q = basix.ufl.quadrature_element("triangle", value_shape=(2,), degree=2)
    S = basix.ufl.element("Lagrange", "triangle", 1, shape=(2, 2), symmetry=True)
    sub_elements = [RT, Q, S] if symmetric else [RT, Q]
    target = ufl.FunctionSpace(
        domain, basix.ufl.mixed_element([DG1, basix.ufl.mixed_element(sub_elements)])
    )
    # The space of each part interpolated into for the reference
    references = {
        "DG1": DG1,
        "RT": RT,
        "Q": Q,
        "S": basix.ufl.element("Lagrange", "triangle", 1, shape=(2, 2)),
    }
    names = ["DG1", "RT", "Q", "S"] if symmetric else ["DG1", "RT", "Q"]
    quadratic = {
        "DG1": ufl.as_vector((x[0] * x[1], x[1] ** 2)),
        "RT": ufl.as_vector((x[0] ** 2, x[1])),
        "Q": ufl.as_vector((x[0] ** 2 * x[1], x[0] * x[1] ** 2)),
        "S": ufl.as_matrix([[x[0] * x[1], x[1] ** 2], [x[1] ** 2, x[0] ** 2]]),
    }
    # In DG1, RT and symmetric P1, and integrated exactly by the quadrature rule
    reproduced = {
        "DG1": ufl.as_vector((1 + 2 * x[0] - x[1], 3 - x[0])),
        "RT": ufl.as_vector((1 + x[0] / 2, -2 + x[1] / 2)),
        "Q": ufl.as_vector((x[0] * x[1], x[1] ** 2)),
        "S": ufl.as_matrix([[1 + x[0], x[1]], [x[1], 2 - x[0]]]),
    }

    def combine(parts):
        """The parts as one vector, their components one after the other."""
        return ufl.as_vector(
            [
                part[i] if part.ufl_shape else part
                for part in (parts[name] for name in names)
                for i in itertools.product(*(range(n) for n in part.ufl_shape))
            ]
        )

    C = ufl.as_vector([float(i + 1) for i in range(combine(quadratic).ufl_shape[0])])
    interpolated = {
        name: ufl.Interpolate(quadratic[name], ufl.FunctionSpace(domain, references[name]))
        for name in names
    }
    forms = [
        ufl.inner(ufl.Interpolate(combine(quadratic), target), C) * ufl.dx,
        ufl.inner(combine(interpolated), C) * ufl.dx,
        ufl.inner(ufl.Interpolate(combine(reproduced), target), C) * ufl.dx,
        ufl.inner(combine(reproduced), C) * ufl.dx,
    ]
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    A = [_tabulate_scalar(module, _integral(module, form), dtype=dtype) for form in compiled]
    for value, value_ref in ((A[0], A[1]), (A[2], A[3])):
        assert abs(value_ref) > 1e-3
        np.testing.assert_allclose(value, value_ref, rtol=500 * np.finfo(dtype).eps)


@pytest.mark.parametrize("kind", ["quadrature", "real"])
def test_interpolate_into_quadrature_or_real_space(compile_args, kind):
    """Interpolating an expression into a vector quadrature or real space.

    The degrees of freedom of a quadrature space are the values at its points,
    so the integral is that of the expression with the same rule. A real space
    is DG0 on a cell: the interpolation is the value at the midpoint of the
    cell, cell by cell, not a value for the whole domain.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    x = ufl.SpatialCoordinate(domain)
    g, c = ufl.as_vector((x[0] * x[1], x[0] * x[1] ** 2)), ufl.as_vector((1.0, 2.0))
    if kind == "quadrature":
        Q = basix.ufl.quadrature_element("triangle", value_shape=(2,), degree=2)
        target = ufl.FunctionSpace(domain, Q)
        forms = [ufl.dot(c, g) * ufl.dx(degree=2)]
    else:
        target = ufl.FunctionSpace(domain, basix.ufl.real_element("triangle", (2,)))
        forms = []
    forms.insert(0, ufl.dot(c, ufl.Interpolate(g, target)) * ufl.dx)
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    A = [_tabulate_scalar(module, _integral(module, form)) for form in compiled]
    if kind == "quadrature":
        A_ref = A[1]
    else:
        # The cell of `_tabulate_scalar`
        vertices = np.array([[0.1, 0.2], [2.0, 0.3], [0.4, 1.7]])
        x0, x1 = vertices.mean(axis=0)
        area = abs(np.linalg.det(vertices[1:] - vertices[0])) / 2
        A_ref = (x0 * x1 + 2 * x0 * x1**2) * area
    assert abs(A_ref) > 1e-3
    np.testing.assert_allclose(A[0], A_ref, rtol=1e-13)
