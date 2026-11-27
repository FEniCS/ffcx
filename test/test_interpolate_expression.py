# Copyright (C) 2026 Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
import sys

import basix
import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.analysis
import ffcx.codegeneration.jit
import ffcx.ir.representation
import ffcx.options
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype

CELL = "triangle"
#: Vertices of an affine triangle that is neither axis aligned nor of unit size
COORDS = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]])
#: The Jacobian of the triangle, with columns v1 - v0 and v2 - v0
JACOBIAN = (COORDS[1:, :2] - COORDS[0, :2]).T
#: Points in the reference triangle to evaluate expressions at
POINTS = np.array([[0.1, 0.1], [0.6, 0.2], [0.25, 0.5], [1 / 3, 1 / 3], [0.05, 0.8]])


def domain():
    return ufl.Mesh(basix.ufl.element("Lagrange", CELL, 1, shape=(2,)))


def space(mesh, family, degree, **kwargs):
    return ufl.FunctionSpace(mesh, basix.ufl.element(family, CELL, degree, **kwargs))


def physical(points):
    """Map reference points to the triangle."""
    return COORDS[0, :2] + points @ JACOBIAN.T


def tolerance(dtype):
    return (
        {"atol": 1e-11, "rtol": 1e-11}
        if np.dtype(dtype).itemsize >= 8
        else {
            "atol": 1e-4,
            "rtol": 1e-4,
        }
    )


def compile_expressions(expressions, compile_args, dtype="float64"):
    compiled, module, _ = ffcx.codegeneration.jit.compile_expressions(
        expressions, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    return compiled, module


def tabulate(
    module,
    expression,
    shape,
    dtype="float64",
    coefficients=None,
    constants=None,
    entity_local_index=None,
):
    """Call the kernel of a compiled expression and reshape what it writes."""
    xdtype = dtype_to_scalar_dtype(dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    A = np.zeros(int(np.prod(shape)), dtype=dtype)
    coords = COORDS.astype(xdtype).flatten()
    w = None if coefficients is None else np.asarray(coefficients, dtype=dtype)
    c = None if constants is None else np.asarray(constants, dtype=dtype)
    entity = None if entity_local_index is None else np.array([entity_local_index], np.int32)
    permutation = np.zeros(1, dtype=np.uint8)

    def pointer(array, ctype):
        return ffi.NULL if array is None else ffi.cast(f"{ctype} *", array.ctypes.data)

    getattr(expression, f"tabulate_tensor_{dtype}")(
        pointer(A, c_type),
        pointer(w, c_type),
        pointer(c, c_type),
        pointer(coords, c_xtype),
        pointer(entity, "int"),
        pointer(permutation, "uint8_t"),
        ffi.NULL,
    )
    return A.reshape(shape)


def lagrange_values(element, dofs, points):
    """Evaluate a function in a (blocked) Lagrange space, ``[point][component]``."""
    block_size = element.block_size
    scalar = element.sub_elements[0] if block_size > 1 else element
    phi = scalar.basix_element.tabulate(0, points)[0][:, :, 0]
    return np.stack([phi @ dofs[c::block_size] for c in range(block_size)], axis=1)


def lagrange_gradient(element, dofs, points):
    """Evaluate the physical gradient of a scalar Lagrange function, ``[point][direction]``."""
    reference = element.basix_element.tabulate(1, points)[1:, :, :, 0]
    return np.einsum("dpi,i->pd", reference, dofs) @ np.linalg.inv(JACOBIAN)


def coefficient_values(element, seed=0):
    return np.random.default_rng(seed).uniform(0.5, 1.5, element.dim)


@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        "float32",
        pytest.param(
            "complex128",
            marks=pytest.mark.skipif(
                sys.platform.startswith("win32"),
                reason="win32 platform does not support C99 _Complex numbers",
            ),
        ),
    ],
)
@pytest.mark.parametrize("target", [("Lagrange", 1, {}), ("Lagrange", 2, {})])
def test_interpolate_coefficient_expression(compile_args, dtype, target):
    """Interpolating an expression of a coefficient, by hand.

    The degrees of freedom of a Lagrange space are point evaluations, so
    ``I(f)(X) = sum_i f(x_i) phi_i(X)`` for the interpolation points ``x_i``.
    """
    mesh = domain()
    V_t = space(mesh, target[0], target[1], **target[2])
    V_u = space(mesh, "Lagrange", 2)
    u = ufl.Coefficient(V_u)
    x = ufl.SpatialCoordinate(mesh)
    expression = ufl.Interpolate(u * u + x[0] * x[1], V_t)

    compiled, module = compile_expressions([(expression, POINTS)], compile_args, dtype)
    (kernel,) = compiled
    assert kernel.num_coefficients == 1
    w = coefficient_values(V_u.ufl_element())
    A = tabulate(module, kernel, (len(POINTS),), dtype, coefficients=w)

    target_element = V_t.ufl_element()
    nodes = target_element.basix_element.points
    u_nodes = lagrange_values(V_u.ufl_element(), w, nodes)[:, 0]
    x_nodes = physical(nodes)
    dofs = u_nodes**2 + x_nodes[:, 0] * x_nodes[:, 1]
    expected = lagrange_values(target_element, dofs, POINTS)[:, 0]
    np.testing.assert_allclose(A, expected, **tolerance(dtype))


@pytest.mark.parametrize(
    "target, expression",
    [
        # Polynomials that lie in the target space are reproduced exactly
        (("Lagrange", 2, {}), lambda x: x[0] ** 2 + x[0] * x[1] - 3 * x[1]),
        (("Lagrange", 1, {"shape": (2,)}), lambda x: ufl.as_vector((x[1] + 1, -2 * x[0]))),
        (
            ("Lagrange", 2, {"shape": (2, 2)}),
            lambda x: ufl.as_matrix([[x[0], x[1] ** 2], [1, x[0] * x[1]]]),
        ),
        # Non-identity interpolation matrices, blocked or Piola mapped
        (
            ("DG", 1, {"shape": (2,), "lagrange_variant": basix.LagrangeVariant.legendre}),
            lambda x: ufl.as_vector((x[0] - x[1], 2 * x[1] + 0.5)),
        ),
        (("N1curl", 1, {}), lambda x: ufl.as_vector((0.3 - x[1], 1.2 + x[0]))),
        # P1 fields plus (x1, -x0) * x0, which is in the second degree space
        (("N1curl", 2, {}), lambda x: ufl.as_vector((x[0] * x[1] + 0.5, x[1] - x[0] ** 2))),
        (("RT", 1, {}), lambda x: ufl.as_vector((0.3 + x[0], -1.2 + x[1]))),
        # Double Piola mapped, which scale with J twice
        (
            ("Regge", 1, {}),
            lambda x: ufl.as_matrix([[1 + x[0], 0.4 - x[1]], [0.4 - x[1], -0.7 + 2 * x[0]]]),
        ),
        (
            ("HHJ", 1, {}),
            lambda x: ufl.as_matrix([[1 + x[0], 0.4 - x[1]], [0.4 - x[1], -0.7 + 2 * x[0]]]),
        ),
    ],
    ids=[
        "P2",
        "P1-vector",
        "P2-tensor",
        "DG1-legendre-vector",
        "N1curl1",
        "N1curl2",
        "RT1",
        "Regge1",
        "HHJ1",
    ],
)
@pytest.mark.parametrize("derivative", [False, True])
def test_interpolate_reproduces_target_space(compile_args, target, expression, derivative):
    """Interpolating a function of the target space gives that function back.

    This holds for any element, so it checks the push forward onto the
    reference cell and the interpolation matrix, blocked or not.
    """
    mesh = domain()
    family, degree, kwargs = target
    V_t = space(mesh, family, degree, **kwargs)
    f = expression(ufl.SpatialCoordinate(mesh))
    interpolated = ufl.Interpolate(f, V_t)
    if derivative:
        f, interpolated = ufl.grad(f), ufl.grad(interpolated)

    compiled, module = compile_expressions([(interpolated, POINTS), (f, POINTS)], compile_args)
    shape = (len(POINTS), *f.ufl_shape)
    A, A_ref = (tabulate(module, kernel, shape) for kernel in compiled)
    assert np.abs(A_ref).max() > 1e-3
    np.testing.assert_allclose(A, A_ref, atol=1e-11)


@pytest.mark.parametrize("family", ["Lagrange", "N1curl"])
def test_interpolate_coefficient_of_target_space(compile_args, family):
    """Interpolating a coefficient into its own space is the identity."""
    mesh = domain()
    V = space(mesh, family, 2)
    u = ufl.Coefficient(V)
    compiled, module = compile_expressions(
        [(ufl.Interpolate(u, V), POINTS), (u, POINTS)], compile_args
    )
    w = coefficient_values(V.ufl_element())
    shape = (len(POINTS), *u.ufl_shape)
    A, A_ref = (tabulate(module, kernel, shape, coefficients=w) for kernel in compiled)
    np.testing.assert_allclose(A, A_ref, atol=1e-12)


def test_interpolate_gradient(compile_args):
    """The gradient of an interpolation is that of the interpolant, by hand."""
    mesh = domain()
    V_t, V_u = space(mesh, "Lagrange", 1), space(mesh, "Lagrange", 2)
    u = ufl.Coefficient(V_u)
    compiled, module = compile_expressions(
        [(ufl.grad(ufl.Interpolate(u * u, V_t)), POINTS)], compile_args
    )
    w = coefficient_values(V_u.ufl_element())
    A = tabulate(module, compiled[0], (len(POINTS), 2), coefficients=w)

    nodes = V_t.ufl_element().basix_element.points
    dofs = lagrange_values(V_u.ufl_element(), w, nodes)[:, 0] ** 2
    expected = lagrange_gradient(V_t.ufl_element(), dofs, POINTS)
    np.testing.assert_allclose(A, expected, atol=1e-12)


@pytest.mark.parametrize(
    "elements",
    [
        (("Lagrange", 2, {}), ("Lagrange", 1, {})),
        (("Lagrange", 1, {}), ("Lagrange", 3, {})),
        (("Lagrange", 2, {"shape": (2,)}), ("Lagrange", 1, {"shape": (2,)})),
        (("N1curl", 2, {}), ("N1curl", 1, {})),
    ],
    ids=["P2-P1", "P1-P3", "vector-P2-P1", "N1curl2-N1curl1"],
)
@pytest.mark.parametrize("derivative", [False, True])
def test_interpolate_argument(compile_args, elements, derivative):
    """Interpolating an argument tabulates the target space against its dofs.

    ``I(phi_j) = sum_i M_ij psi_i`` for the reference interpolation operator
    ``M``, so the values must be those of the target space's argument times
    ``M``.
    """
    mesh = domain()
    (fam_a, deg_a, kw_a), (fam_t, deg_t, kw_t) = elements
    V_a, V_t = space(mesh, fam_a, deg_a, **kw_a), space(mesh, fam_t, deg_t, **kw_t)
    interpolated, reference = ufl.Interpolate(ufl.TestFunction(V_a), V_t), ufl.TestFunction(V_t)
    if derivative:
        interpolated, reference = ufl.grad(interpolated), ufl.grad(reference)

    compiled, module = compile_expressions(
        [(interpolated, POINTS), (reference, POINTS)], compile_args
    )
    element_a, element_t = V_a.ufl_element(), V_t.ufl_element()
    assert compiled[0].rank == 1
    components = int(np.prod(interpolated.ufl_shape))
    A = tabulate(module, compiled[0], (len(POINTS), components, element_a.dim))
    A_ref = tabulate(module, compiled[1], (len(POINTS), components, element_t.dim))

    def base(element):
        return element.sub_elements[0] if element.block_size > 1 else element

    operator = basix.compute_interpolation_operator(
        base(element_a).basix_element, base(element_t).basix_element
    )
    operator = np.kron(operator, np.eye(element_a.block_size))
    assert np.abs(A).max() > 1e-3
    np.testing.assert_allclose(A, A_ref @ operator, atol=1e-12)


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("target_degree", [1, 2])
def test_interpolate_linear_in_argument(compile_args, dtype, target_degree):
    """An expression linear in an argument needs its table built for each cell.

    ``A[q, j] = sum_i (u phi_j)(x_i) psi_i(X_q)`` for the target's nodes ``x_i``.
    """
    mesh = domain()
    V_t, V_a = space(mesh, "Lagrange", target_degree), space(mesh, "Lagrange", 2)
    u, v = ufl.Coefficient(V_a), ufl.TestFunction(V_a)
    x = ufl.SpatialCoordinate(mesh)
    compiled, module = compile_expressions(
        [(ufl.Interpolate((u + x[0]) * v, V_t), POINTS)], compile_args, dtype
    )
    element_a, element_t = V_a.ufl_element(), V_t.ufl_element()
    w = coefficient_values(element_a)
    A = tabulate(module, compiled[0], (len(POINTS), element_a.dim), dtype, coefficients=w)

    nodes = element_t.basix_element.points
    phi = element_a.basix_element.tabulate(0, nodes)[0][:, :, 0]
    scale = lagrange_values(element_a, w, nodes)[:, 0] + physical(nodes)[:, 0]
    psi = element_t.basix_element.tabulate(0, POINTS)[0][:, :, 0]
    expected = psi @ (scale[:, None] * phi)
    np.testing.assert_allclose(A, expected, **tolerance(dtype))


def test_interpolate_linear_in_argument_blocked(compile_args):
    """A table built for each cell, for a blocked target and argument."""
    mesh = domain()
    V = space(mesh, "Lagrange", 1, shape=(2,))
    u, v = ufl.Coefficient(V), ufl.TestFunction(V)
    expression = ufl.as_vector((u[0] * v[0], u[1] * v[1] + u[0] * v[1]))
    compiled, module = compile_expressions([(ufl.Interpolate(expression, V), POINTS)], compile_args)
    element = V.ufl_element()
    w = coefficient_values(element)
    A = tabulate(module, compiled[0], (len(POINTS), 2, element.dim), coefficients=w)

    # The target is P1, so the interpolant is linear between its vertex values
    nodes = element.sub_elements[0].basix_element.points
    phi = element.sub_elements[0].basix_element.tabulate(0, nodes)[0][:, :, 0]
    u_nodes = lagrange_values(element, w, nodes)
    values = np.zeros((len(nodes), 2, element.dim))  # [node][component][dof]
    for j in range(element.dim):
        node, component = divmod(j, 2)
        if component == 0:
            values[:, 0, j] = u_nodes[:, 0] * phi[:, node]
        else:
            values[:, 1, j] = (u_nodes[:, 1] + u_nodes[:, 0]) * phi[:, node]
    psi = element.sub_elements[0].basix_element.tabulate(0, POINTS)[0][:, :, 0]
    np.testing.assert_allclose(A, np.einsum("qn,ncj->qcj", psi, values), atol=1e-12)


def finite_difference(evaluate, w, step=1e-6):
    """Central differences of ``evaluate`` with respect to each entry of ``w``."""
    columns = []
    for k in range(len(w)):
        direction = np.zeros_like(w)
        direction[k] = step
        columns.append((evaluate(w + direction) - evaluate(w - direction)) / (2 * step))
    return np.stack(columns, axis=-1)


@pytest.mark.parametrize("ngrads", [0, 1])
@pytest.mark.parametrize("block_size", [1, 2])
def test_derivative_of_interpolation(compile_args, ngrads, block_size):
    """Differentiating an expression through an interpolation.

    ``d/du I(g(u))[v] = I(dg/du[v])``, which is linear in the argument and so
    needs its table built for each cell. Checked against finite differences of
    the expression itself.

    UFL takes an interpolation that is the whole expression for a form, whose
    coargument is numbered 0, so the derivative is taken in the direction of an
    argument numbered 1.
    """
    mesh = domain()
    shape = () if block_size == 1 else (block_size,)
    V_a = space(mesh, "Lagrange", 1, shape=shape)
    V_t = space(mesh, "Lagrange", 2 if block_size == 1 else 1, shape=shape)
    u, v = ufl.Coefficient(V_a), ufl.TrialFunction(V_a)
    x = ufl.SpatialCoordinate(mesh)
    if block_size == 1:
        g = u * u * x[0]
    else:
        g = ufl.as_vector((u[0] * u[1], u[1] * u[1] + x[1]))
    expression = ufl.Interpolate(g, V_t)
    if ngrads:
        expression = ufl.grad(expression)
    derivative = ufl.derivative(expression, u, v)

    compiled, module = compile_expressions(
        [(expression, POINTS), (derivative, POINTS)], compile_args
    )
    n = V_a.ufl_element().dim
    shape = (len(POINTS), int(np.prod(expression.ufl_shape)))
    w = coefficient_values(V_a.ufl_element())
    A = tabulate(module, compiled[1], (*shape, n), coefficients=w)
    expected = finite_difference(lambda w: tabulate(module, compiled[0], shape, coefficients=w), w)
    assert np.abs(expected).max() > 1e-3
    np.testing.assert_allclose(A, expected, atol=1e-7)


@pytest.mark.parametrize("ngrads", [0, 1])
def test_second_derivative_of_interpolation(compile_args, ngrads):
    """Two derivatives give an interpolation bilinear in two arguments.

    Its table is a dense block, with one axis per argument, built for each cell.
    Checked against finite differences of the first derivative.

    The interpolation is scaled so that it is not the whole expression, which
    UFL would take for a form whose coargument clashes with an argument.
    """
    mesh = domain()
    V_a, V_t = space(mesh, "Lagrange", 1), space(mesh, "Lagrange", 2)
    u = ufl.Coefficient(V_a)
    v0, v1 = ufl.TestFunction(V_a), ufl.TrialFunction(V_a)
    x = ufl.SpatialCoordinate(mesh)
    expression = ufl.Interpolate(u**3, V_t)
    if ngrads:
        expression = ufl.grad(expression)[0]
    expression = expression * (1 + x[0])
    first = ufl.derivative(expression, u, v0)
    second = ufl.derivative(first, u, v1)

    compiled, module = compile_expressions([(first, POINTS), (second, POINTS)], compile_args)
    assert compiled[1].rank == 2
    n = V_a.ufl_element().dim
    w = coefficient_values(V_a.ufl_element())
    A = tabulate(module, compiled[1], (len(POINTS), 1, n, n), coefficients=w)
    expected = finite_difference(
        lambda w: tabulate(module, compiled[0], (len(POINTS), 1, n), coefficients=w), w
    )
    assert np.abs(expected).max() > 1e-3
    np.testing.assert_allclose(A, expected, atol=1e-7)


def test_interpolation_coefficient_order(compile_args):
    """Coefficients inside and outside of an interpolation are passed once each.

    The kernel takes the coefficients of the interpolated expression from the
    ones passed to the expression, in the order the expression gives them.
    """
    mesh = domain()
    V1, V2 = space(mesh, "Lagrange", 1), space(mesh, "Lagrange", 2)
    a = ufl.Coefficient(V2)  # only inside the interpolation
    b = ufl.Coefficient(V1)  # only outside
    c = ufl.Coefficient(V2)  # both
    expression = b * c * ufl.Interpolate(a * c, V2)

    compiled, module = compile_expressions([(expression, POINTS)], compile_args)
    (kernel,) = compiled
    assert kernel.num_coefficients == 3
    positions = [kernel.original_coefficient_positions[i] for i in range(3)]
    assert positions == [0, 1, 2]  # a, b, c, as extract_coefficients gives them

    values = {f: coefficient_values(f.ufl_element(), seed) for seed, f in enumerate((a, b, c))}
    A = tabulate(module, kernel, (len(POINTS),), coefficients=np.concatenate(list(values.values())))
    nodes = V2.ufl_element().basix_element.points

    def at(f, points):
        return lagrange_values(f.ufl_element(), values[f], points)[:, 0]

    interpolated = lagrange_values(V2.ufl_element(), at(a, nodes) * at(c, nodes), POINTS)[:, 0]
    np.testing.assert_allclose(A, at(b, POINTS) * at(c, POINTS) * interpolated, atol=1e-12)


def test_interpolation_constants(compile_args):
    """The constants of an interpolated expression are passed to its kernel.

    They are passed in the order the kernel reads them, which is not that of
    the expression holding the interpolation.
    """
    mesh = domain()
    V = space(mesh, "Lagrange", 2)
    u = ufl.Coefficient(V)
    c0, c1, c2 = ufl.Constant(mesh), ufl.Constant(mesh, shape=(2,)), ufl.Constant(mesh)
    expression = c0 * ufl.Interpolate(c2 * u + c1[1], V) + c1[0]

    compiled, module = compile_expressions([(expression, POINTS)], compile_args)
    (kernel,) = compiled
    assert kernel.num_constants == 3
    constants = np.array([2.0, 3.0, -0.5, 0.25])  # c0, c1[0], c1[1], c2
    w = coefficient_values(V.ufl_element())
    A = tabulate(module, kernel, (len(POINTS),), coefficients=w, constants=constants)
    # The target is the space of u, so the interpolation reproduces it
    u_values = lagrange_values(V.ufl_element(), w, POINTS)[:, 0]
    np.testing.assert_allclose(A, 2.0 * (0.25 * u_values - 0.5) + 3.0, atol=1e-12)


def test_interpolations_sharing_a_kernel(compile_args):
    """Interpolations that differ only in their coefficients share a kernel.

    Each is still evaluated with its own coefficients, and the shared kernel is
    generated once, as is the kernel of an interpolation in two expressions.
    """
    mesh = domain()
    V1, V2 = space(mesh, "Lagrange", 1), space(mesh, "Lagrange", 2)
    u, w = ufl.Coefficient(V2), ufl.Coefficient(V2)
    first = ufl.Interpolate(u * u, V1) - 2 * ufl.Interpolate(w * w, V1)
    second = ufl.grad(ufl.Interpolate(u * u, V1))[0]

    analysis = ffcx.analysis.analyze_ufl_objects([(first, POINTS), (second, POINTS)], np.float64)
    options = ffcx.options.get_options({"scalar_type": "float64"})
    ir = ffcx.ir.representation.compute_ir(analysis, {}, "", options, None)
    assert len(analysis.sub_kernels) == 3  # two proxies in each expression
    assert len(ir.sub_kernels) == 1  # one kernel, shared
    assert len(ir.expressions) == 2

    compiled, module = compile_expressions([(first, POINTS), (second, POINTS)], compile_args)
    u_values = coefficient_values(V2.ufl_element(), 0)
    w_values = coefficient_values(V2.ufl_element(), 1)
    A = tabulate(module, compiled[0], (len(POINTS),), coefficients=np.r_[u_values, w_values])
    nodes = V1.ufl_element().basix_element.points

    def interpolant(values):
        """The dofs of the interpolation of the square of a function in V2."""
        return lagrange_values(V2.ufl_element(), values, nodes)[:, 0] ** 2

    expected = lagrange_values(
        V1.ufl_element(), interpolant(u_values) - 2 * interpolant(w_values), POINTS
    )[:, 0]
    np.testing.assert_allclose(A, expected, atol=1e-12)

    B = tabulate(module, compiled[1], (len(POINTS),), coefficients=u_values)
    gradient = lagrange_gradient(V1.ufl_element(), interpolant(u_values), POINTS)
    np.testing.assert_allclose(B, gradient[:, 0], atol=1e-12)


def test_nested_interpolation(compile_args):
    """An interpolation of an expression holding an interpolation, by hand."""
    mesh = domain()
    V1, V2 = space(mesh, "Lagrange", 1), space(mesh, "Lagrange", 2)
    u = ufl.Coefficient(V2)
    expression = ufl.Interpolate(ufl.Interpolate(u * u, V1) * u, V2)

    compiled, module = compile_expressions([(expression, POINTS)], compile_args)
    assert compiled[0].num_coefficients == 1
    w = coefficient_values(V2.ufl_element())
    A = tabulate(module, compiled[0], (len(POINTS),), coefficients=w)

    def u_at(points):
        return lagrange_values(V2.ufl_element(), w, points)[:, 0]

    nodes1, nodes2 = (V.ufl_element().basix_element.points for V in (V1, V2))
    inner = u_at(nodes1) ** 2
    outer = lagrange_values(V1.ufl_element(), inner, nodes2)[:, 0] * u_at(nodes2)
    expected = lagrange_values(V2.ufl_element(), outer, POINTS)[:, 0]
    np.testing.assert_allclose(A, expected, atol=1e-12)


@pytest.mark.parametrize("facet", [0, 1, 2])
def test_interpolation_on_facet(compile_args, facet):
    """An interpolation in an expression evaluated on a facet.

    The interpolated expression is evaluated in the whole cell, so the values
    on a facet must be those in the cell at the same points.
    """
    mesh = domain()
    V1, V2 = space(mesh, "Lagrange", 1), space(mesh, "Lagrange", 2)
    u = ufl.Coefficient(V2)
    expression = ufl.Interpolate(u * u, V1) * ufl.Interpolate(ufl.TestFunction(V2), V1)

    facet_points = np.array([[0.1], [0.45], [0.8]])
    triangle = basix.CellType.triangle
    vertices = basix.geometry(triangle)[basix.topology(triangle)[1][facet]]
    cell_points = vertices[0] + facet_points * (vertices[1] - vertices[0])

    compiled, module = compile_expressions(
        [(expression, facet_points), (expression, cell_points)], compile_args
    )
    w = coefficient_values(V2.ufl_element())
    shape = (len(facet_points), V2.ufl_element().dim)
    A = tabulate(module, compiled[0], shape, coefficients=w, entity_local_index=facet)
    A_ref = tabulate(module, compiled[1], shape, coefficients=w)
    np.testing.assert_allclose(A, A_ref, atol=1e-12)


def test_unsupported_interpolations():
    """Interpolations that cannot be compiled are reported as such."""
    mesh = domain()
    V = space(mesh, "Lagrange", 1)
    u, v = ufl.Coefficient(V), ufl.TestFunction(V)

    # A table built for each cell is only contracted for the cell itself
    with pytest.raises(NotImplementedError, match="only linear in the argument"):
        ffcx.codegeneration.jit.compile_expressions(
            [(ufl.Interpolate(u * v, V), np.array([[0.5]]))]
        )

    # Mapping an expression that is only linear in the argument onto the
    # reference cell of a Piola mapped space is not implemented
    W = space(mesh, "N1curl", 1)
    with pytest.raises(NotImplementedError, match="space is only supported"):
        ffcx.codegeneration.jit.compile_expressions(
            [(ufl.Interpolate(u * ufl.TestFunction(W), W), POINTS)]
        )
    # Including one whose argument is shaped differently from the space
    with pytest.raises(NotImplementedError, match="space is only supported"):
        ffcx.codegeneration.jit.compile_expressions([(ufl.Interpolate(ufl.grad(v), W), POINTS)])

    # The dofs of a proxy coefficient are built from a single basix element
    x = ufl.SpatialCoordinate(mesh)
    symmetric = space(mesh, "Lagrange", 1, shape=(2, 2), symmetry=True)
    mixed = ufl.FunctionSpace(
        mesh,
        basix.ufl.mixed_element(
            [basix.ufl.element("Lagrange", CELL, 1), basix.ufl.element("N1curl", CELL, 1)]
        ),
    )
    for target, expression in [
        (symmetric, ufl.as_matrix([[x[0], x[1]], [x[1], 1.0]])),
        (mixed, ufl.as_vector((x[0], x[1], 1.0))),
    ]:
        with pytest.raises(NotImplementedError, match="space is not supported"):
            ffcx.codegeneration.jit.compile_expressions(
                [(ufl.Interpolate(expression, target), POINTS)]
            )

    # As is the kernel evaluating an expression that is only linear in the
    # argument at the interpolation points of the target
    quadrature = ufl.FunctionSpace(mesh, basix.ufl.quadrature_element(CELL, degree=2))
    real = ufl.FunctionSpace(mesh, basix.ufl.real_element(CELL, ()))
    for target in [quadrature, real]:
        with pytest.raises(NotImplementedError, match="space is not supported"):
            ffcx.codegeneration.jit.compile_forms([ufl.Interpolate(u * v, target) * ufl.dx])


def test_sub_kernel_targets():
    """What needs a sub-kernel of its own, and what it reads."""
    mesh = domain()
    V1, V2 = space(mesh, "Lagrange", 1), space(mesh, "Lagrange", 2)
    u, z = ufl.Coefficient(V2), ufl.Coefficient(V1)
    c = ufl.Constant(mesh)
    v = ufl.TestFunction(V2)
    expression = (
        z
        * ufl.Interpolate(c * u * u, V1)  # a proxy coefficient
        * ufl.Interpolate(v, V1)  # a table known at compile time
        + ufl.Interpolate(u * v, V1)  # a table built for each cell
    )
    analysis = ffcx.analysis.analyze_ufl_objects([(expression, POINTS)], np.float64)
    processed, _, _ = analysis.expressions[-1]

    interpolations = ffcx.analysis.sub_kernel_targets([processed])
    assert len(interpolations) == 2
    proxy, table = interpolations
    assert isinstance(proxy, ffcx.analysis.ProxyCoefficient)
    assert isinstance(table, ufl.Interpolate)
    assert proxy.is_interpolation
    assert {sub_kernel.target for sub_kernel in analysis.sub_kernels} == {proxy, table}
    for sub_kernel in analysis.sub_kernels:
        assert sub_kernel.expression is not None and sub_kernel.form_data is None

    assert ffcx.analysis.sub_kernel_coefficients(proxy) == [u]
    assert ffcx.analysis.sub_kernel_constants(proxy) == [c]
    assert ffcx.analysis.sub_kernel_coefficients(table) == [u]
    # The proxy coefficient stands for what it is evaluated from
    assert ffcx.analysis.required_coefficients(processed) == [u, z]
    assert ffcx.analysis.required_constants(processed) == [c]
