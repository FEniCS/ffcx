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


@pytest.mark.parametrize("degrees", [(2, 1), (1, 2), (3, 2)])
def test_interpolate_argument_tensor_product(compile_args, degrees):
    """Interpolating an argument between tensor product spaces on a quadrilateral.

    Sum factorisation is on by default for these elements. The interpolated
    table includes the interpolation operator, which the 1D tables of the
    target element leave out, so it must not be sum factorised.
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
    J = ufl.Interpolate(ufl.TrialFunction(V_a), V_t) * ufl.dx
    J_ref = ufl.TrialFunction(V_t) * ufl.dx

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    coords = np.array(
        [[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0], [2.3, 1.9, 0.0]], dtype=np.float64
    ).flatten()
    ffi = module.ffi

    def tabulate(form, size):
        offsets = form.form_integral_offsets
        integral = form.form_integrals[offsets[module.lib.cell]]
        A = np.zeros(size, dtype=dtype)
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

    A = tabulate(compiled_forms[0], element_a.dim)
    A_ref = tabulate(compiled_forms[1], element_t.dim)
    interpolation = basix.compute_interpolation_operator(
        element_a.basix_element, element_t.basix_element
    )
    np.testing.assert_allclose(A, interpolation.T @ A_ref, atol=1e-13)


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


@pytest.mark.parametrize("do_cancel_jacobian_products", [False, True])
@pytest.mark.parametrize("family", ["Lagrange", "N1curl"])
def test_interpolate_argument_is_compile_time(family, do_cancel_jacobian_products):
    """Interpolating an argument itself must not need a table built per cell.

    The operand is compared against the argument as it stands after the pull
    backs, which lower the geometry of a Piola map, and after cancelling
    Jacobian products, which removes its component tensors, so the comparison
    has to be made up to both.
    """
    cell = "triangle"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    V_a = ufl.FunctionSpace(domain, basix.ufl.element(family, cell, 2))
    V_t = ufl.FunctionSpace(domain, basix.ufl.element(family, cell, 1))

    interpolated = ufl.Interpolate(ufl.TrialFunction(V_a), V_t)
    integrand = interpolated
    if integrand.ufl_shape:
        integrand = sum(integrand[i] for i in range(integrand.ufl_shape[0]))

    form = ufl.algorithms.compute_form_data(
        integrand * ufl.dx,
        do_apply_function_pullbacks=True,
        do_apply_geometry_lowering=True,
        preserve_geometry_types=(ufl.classes.Jacobian,),
        do_cancel_jacobian_products=do_cancel_jacobian_products,
    )
    interpolations = [
        node
        for integral_data in form.integral_data
        for integral in integral_data.integrals
        for node in ufl.corealg.traversal.unique_pre_traversal(integral.integrand())
        if isinstance(node, ufl.Interpolate)
    ]
    assert interpolations
    for node in interpolations:
        assert ffcx.analysis.interpolated_argument(node) is not None
        assert not ffcx.analysis.interpolation_has_runtime_table(node)


def _integral(module, form, subdomain_id=-1):
    """Get the cell integral of a compiled form over a subdomain."""
    cell = module.lib.cell
    offsets = form.form_integral_offsets
    for i in range(offsets[cell], offsets[cell + 1]):
        if form.form_integral_ids[i] == subdomain_id:
            return form.form_integrals[i]
    raise KeyError(subdomain_id)


def _tabulate_scalar(module, integral, coefficients=None, constants=None):
    """Tabulate a functional on an affine triangle."""
    ffi = module.ffi
    coords = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]]).flatten()
    A = np.zeros(1)

    def pointer(array):
        return ffi.NULL if array is None else ffi.cast("double *", array.ctypes.data)

    integral.tabulate_tensor_float64(
        pointer(A),
        pointer(coefficients),
        pointer(constants),
        pointer(coords),
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


def test_interpolate_into_blocked_space_with_interpolation_matrix(compile_args):
    """A blocked target space whose interpolation matrix is not the identity.

    Each block component is interpolated with the matrix of the scalar
    sub-element. A linear field lies in the space, so it is reproduced.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    element = basix.ufl.element(
        "DG", "triangle", 1, shape=(2,), lagrange_variant=basix.LagrangeVariant.legendre
    )
    assert not element.basix_element.interpolation_is_identity
    V = ufl.FunctionSpace(domain, element)
    x = ufl.SpatialCoordinate(domain)
    f = ufl.as_vector((x[0] - 2 * x[1], 0.5 + x[1]))
    forms = [
        ufl.dot(ufl.Interpolate(f, V), ufl.as_vector((1.0, 3.0))) * ufl.dx,
        ufl.dot(f, ufl.as_vector((1.0, 3.0))) * ufl.dx,
    ]
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    A, A_ref = (_tabulate_scalar(module, _integral(module, form)) for form in compiled)
    assert abs(A_ref) > 1e-3
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
    identity, so each form must match the one without the interpolation. The
    dofs of a proxy coefficient are only built for the '+' cell, so '-' is
    reported as not implemented.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    v, u = ufl.TestFunction(V), ufl.Coefficient(V)
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

    def tabulate(form, size):
        (compiled,), module, _ = ffcx.codegeneration.jit.compile_forms(
            [form], cffi_extra_compile_args=compile_args
        )
        ffi = module.ffi
        offsets = compiled.form_integral_offsets
        integral = compiled.form_integrals[offsets[module.lib.interior_facet]]
        A = np.zeros(size)
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

    A = tabulate(ufl.Interpolate(u, V)("+") * (1 + x[0]) * ufl.dS, 1)
    A_ref = tabulate(u("+") * (1 + x[0]) * ufl.dS, 1)
    np.testing.assert_allclose(A, A_ref, rtol=1e-13)

    with pytest.raises(NotImplementedError, match="restricted to '-'"):
        tabulate(ufl.Interpolate(u, V)("-") * ufl.dS, 1)


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
    (proxies,) = (
        integral.interpolation.proxies
        for integral in ir.integrals
        if integral.interpolation.proxies
    )
    assert len(proxies) == 1
