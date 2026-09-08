import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
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
        [J, J_ref],
        options={"scalar_type": dtype},
        # cache_dir=".ffcx_cache",
        cffi_extra_compile_args=[],  # compile_args,
        visualise=False,
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
def test_interpolate_derivative_in_form(compile_args, dtype, ngrads):
    """Differentiating a form holding an interpolation applies the chain rule.

    ``d/du F(I(u))[v] = <dF/dI, I(v)>``, and interpolation is linear, so the
    tensor must be the reference interpolation operator applied to the tensor of
    the same form posed on the target space.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element_a = basix.ufl.element("Lagrange", cell, 1)
    element_t = basix.ufl.element("Lagrange", cell, 2)
    V_a = ufl.FunctionSpace(domain, element_a)
    V_t = ufl.FunctionSpace(domain, element_t)

    u, du = ufl.Coefficient(V_a), ufl.TestFunction(V_a)
    v_t = ufl.TrialFunction(V_t)

    def differentiate(expression):
        """Take ``ngrads`` gradients and index back down to a scalar."""
        for _ in range(ngrads):
            expression = ufl.grad(expression)
        return expression[(0,) * ngrads] if ngrads else expression

    J = ufl.derivative(differentiate(ufl.Interpolate(u, V_t)) * ufl.dx, u, du)
    J_ref = differentiate(v_t) * ufl.dx

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

    A = tabulate(compiled_forms[0], element_a.dim)
    A_ref = tabulate(compiled_forms[1], element_t.dim)

    interpolation = basix.compute_interpolation_operator(
        element_a.basix_element, element_t.basix_element
    )
    tol = np.finfo(dtype).eps * 100
    assert np.abs(A).max() > tol
    np.testing.assert_allclose(A, interpolation.T @ A_ref, atol=tol)


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
