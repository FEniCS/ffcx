# Copyright (C) 2026 Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
import itertools

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.analysis
import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


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


@pytest.mark.parametrize("do_cancel_jacobian_products", [False, True])
@pytest.mark.parametrize("scalar_type", [np.float64, np.complex128])
@pytest.mark.parametrize("measure", ["dx", "ds", "expression"])
@pytest.mark.parametrize("family", ["Lagrange", "N1curl", "RT", "Regge", "HHJ"])
def test_interpolate_argument_is_compile_time(
    family, measure, scalar_type, do_cancel_jacobian_products
):
    """Interpolating an argument itself must not need a table built per cell.

    UFL leaves the argument mapped by its own pull back and the inverse of the
    target's, with the geometry of a Piola map lowered, and reduces it to its
    reference value only if Jacobian products are cancelled. The analysis
    recognises it either way, and rewrites it into its reference value.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    V_a = ufl.FunctionSpace(domain, basix.ufl.element(family, cell, 2))
    V_t = ufl.FunctionSpace(domain, basix.ufl.element(family, cell, 1))
    argument = ufl.TestFunction(V_a)

    interpolated = ufl.Interpolate(argument, V_t)
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


def test_interpolation_on_interior_facet(compile_args):
    """Interpolations restricted to either side of an interior facet.

    Interpolating into the space the argument lives in is the identity, so each
    form must match the one without the interpolation. A bilinear form
    restricts its two interpolations to each combination of sides.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    v, u = ufl.TestFunction(V), ufl.TrialFunction(V)
    x = ufl.SpatialCoordinate(domain)

    # Two triangles sharing the facet (1, 0) -> (0, 1), local facet 0 of both,
    # with the same vertex order so neither side is permuted.
    coords = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float64
    ).flatten()
    entity = np.array([0, 0], dtype=np.int32)
    perm = np.array([0, 0], dtype=np.uint8)

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
            ffi.NULL,
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
