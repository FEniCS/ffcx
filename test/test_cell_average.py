# Copyright (C) 2026 Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""The cell average of an expression, evaluated by sub-kernels.

``cell_avg(f)`` is the ratio of two integrals over the cell, of ``f`` and of
one, each assembled by the cell kernel of a linear form tested with the basis
function of a P0 space, which is one. The kernel that holds the average calls
the two sub-kernels and divides.
"""

import basix
import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.analysis
import ffcx.codegeneration.jit
import ffcx.compiler
import ffcx.ir.representation
import ffcx.options
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


def curved_nodes(cell, degree):
    """Nodes of a curved cell: the reference nodes of its coordinate element, mapped nonlinearly."""
    X = basix.ufl.element("Lagrange", cell, degree).basix_element.points
    x = np.column_stack(
        [
            2.0 * X[:, 0] + 0.5 * X[:, 1] ** 3,
            2.0 * X[:, 1] + 0.5 * X[:, 0] ** 2 * X[:, 1],
            np.zeros(len(X)),
        ]
    )
    return x


#: Cells to average on, each as the cell type, the degree of its coordinate
#: element and the coordinates of its nodes: a triangle that is neither axis
#: aligned nor of unit size, a quadrilateral that is not affine, and triangles
#: whose edges are curved.
GEOMETRIES = {
    "triangle": ("triangle", 1, np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]])),
    "quadrilateral": (
        "quadrilateral",
        1,
        np.array([[0.0, 0.0, 0.0], [2.0, 0.2, 0.0], [0.1, 1.5, 0.0], [2.5, 2.0, 0.0]]),
    ),
    "P2 triangle": ("triangle", 2, curved_nodes("triangle", 2)),
    "P3 triangle": ("triangle", 3, curved_nodes("triangle", 3)),
}


def tabulate(
    module,
    compiled_form,
    shape,
    coefficients=(),
    constants=(),
    coords=GEOMETRIES["triangle"][2],
    dtype="float64",
    integral_type="cell",
    entity=None,
    permutation=None,
):
    """Assemble the one integral of a compiled form, of the given type, on one cell."""
    offsets = compiled_form.form_integral_offsets
    kind = getattr(module.lib, integral_type)
    assert offsets[kind + 1] - offsets[kind] == 1
    integral = compiled_form.form_integrals[offsets[kind]]

    xdtype = dtype_to_scalar_dtype(dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    A = np.zeros(shape, dtype=dtype)
    w = np.asarray(coefficients, dtype=dtype)
    c = np.asarray(constants, dtype=dtype)
    x = np.asarray(coords, dtype=xdtype).flatten()
    entity = ffi.NULL if entity is None else ffi.cast("int *", np.asarray(entity).ctypes.data)
    permutation = (
        ffi.NULL
        if permutation is None
        else ffi.cast("uint8_t *", np.asarray(permutation).ctypes.data)
    )
    getattr(integral, f"tabulate_tensor_{dtype}")(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", x.ctypes.data),
        entity,
        permutation,
        ffi.NULL,
    )
    return A


def compile_forms(forms, compile_args, dtype="float64"):
    """Compile forms into one module."""
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    return compiled, module


def spaces(cell="triangle", degree=1):
    """A mesh, with a P2 space for the coefficient and a P1 space for the test function."""
    mesh = ufl.Mesh(basix.ufl.element("Lagrange", cell, degree, shape=(2,)))
    V = ufl.FunctionSpace(mesh, basix.ufl.element("Lagrange", cell, 2))
    Q = ufl.FunctionSpace(mesh, basix.ufl.element("Lagrange", cell, 1))
    return mesh, V, Q


def integral_forms(f, mesh):
    """The forms whose vectors are the integrals of an expression and of one over the cell.

    The expression is tested with the basis function of a P0 space, which is
    one, as the sub-kernels of a cell average are built.
    """
    P0 = ufl.FunctionSpace(mesh, basix.ufl.element("DG", mesh.ufl_cell().cellname, 0))
    q = ufl.TestFunction(P0)
    if f.ufl_shape:
        P0_f = ufl.FunctionSpace(
            mesh, basix.ufl.element("DG", mesh.ufl_cell().cellname, 0, shape=f.ufl_shape)
        )
        numerator = ufl.inner(f, ufl.TestFunction(P0_f)) * ufl.dx
    else:
        numerator = f * q * ufl.dx
    return numerator, q * ufl.dx


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("geometry", list(GEOMETRIES))
@pytest.mark.parametrize("vector", [False, True])
def test_cell_average(compile_args, dtype, geometry, vector):
    """The cell average of an expression of a coefficient, against quadrature by hand."""
    cell, degree, coords = GEOMETRIES[geometry]
    mesh, V, Q = spaces(cell, degree)
    u = ufl.Coefficient(V)
    v = ufl.TestFunction(Q)
    x = ufl.SpatialCoordinate(mesh)
    if vector:
        # Reads the gradient through the Jacobian in the sub-kernel
        form = ufl.inner(ufl.cell_avg(u * ufl.grad(u)), ufl.grad(v)) * ufl.dx
    else:
        form = ufl.cell_avg(u * u + x[0] * x[1]) * v * ufl.dx

    dofs = np.random.default_rng(0).uniform(0.5, 1.5, V.ufl_element().dim)
    (compiled,), module = compile_forms([form], compile_args, dtype)
    A = tabulate(module, compiled, Q.ufl_element().dim, dofs, coords=coords, dtype=dtype)

    # The geometry at the points of a quadrature rule exact for the integrands
    points, weights = basix.make_quadrature(basix.CellType[cell], 12)
    chi = basix.ufl.element("Lagrange", cell, degree).basix_element.tabulate(1, points)
    x_q = chi[0, :, :, 0] @ coords[:, :2]
    J_q = np.einsum("dqi,ix->qxd", chi[1:, :, :, 0], coords[:, :2])
    detJ_q = np.abs(np.linalg.det(J_q))
    invJ_q = np.linalg.inv(J_q)

    def gradient(reference):
        """Push forward reference gradients ``[direction][point][function]``."""
        return np.einsum("dqi,qdx->qix", reference, invJ_q)

    phi = V.ufl_element().basix_element.tabulate(1, points)
    u_q = phi[0, :, :, 0] @ dofs
    if vector:
        f_q = u_q[:, None] * np.einsum("qix,i->qx", gradient(phi[1:, :, :, 0]), dofs)
    else:
        f_q = u_q**2 + x_q[:, 0] * x_q[:, 1]
    average = (weights * detJ_q) @ f_q / (weights @ detJ_q)

    psi = Q.ufl_element().basix_element.tabulate(1, points)
    if vector:
        expected = np.einsum("q,qix,x->i", weights * detJ_q, gradient(psi[1:, :, :, 0]), average)
    else:
        expected = average * ((weights * detJ_q) @ psi[0, :, :, 0])

    tol = 1e-4 if np.dtype(dtype).itemsize < 8 else 1e-11
    np.testing.assert_allclose(A, expected, rtol=tol, atol=tol)


@pytest.mark.parametrize("shape", [(), (2,), (2, 2)])
def test_cell_average_is_ratio_of_integrals(compile_args, shape):
    """A cell average is the integral of the expression divided by the volume.

    Each is the vector of a linear form tested with the basis function of a
    P0 space, and the average is read as a P0 function: the form that holds it
    must match the forms of the integrals assembled on their own, to rounding.
    The expression reads a constant, packed for the sub-kernel, and is a
    tensor of the given shape.
    """
    mesh, V, Q = spaces()
    u = ufl.Coefficient(V)
    c = ufl.Constant(mesh, shape=(2,))
    v = ufl.TestFunction(Q)
    x = ufl.SpatialCoordinate(mesh)
    if shape == ():
        f = c[0] * u * u + c[1] * x[1]
    elif shape == (2,):
        f = c[0] * u * ufl.grad(u) + c[1] * x
    else:
        f = c[0] * ufl.outer(ufl.grad(u), x) + c[1] * ufl.sym(ufl.outer(x, ufl.grad(u)))
    average = ufl.cell_avg(f)
    if shape == ():
        form = average * v * ufl.dx
    else:
        weight = ufl.grad(v) if shape == (2,) else ufl.outer(ufl.grad(v), x)
        form = ufl.inner(average, weight) * ufl.dx
    numerator, volume = integral_forms(f, mesh)
    # The average is a P0 function, so this form with it as a coefficient is
    # the form that the cell average is read through
    P0 = ufl.FunctionSpace(
        mesh, basix.ufl.element("DG", "triangle", 0, shape=shape if shape else None)
    )
    p = ufl.Coefficient(P0)
    if shape == ():
        weighted = p * v * ufl.dx
    else:
        weighted = ufl.inner(p, weight) * ufl.dx

    forms = [form, numerator, volume, weighted]
    compiled, module = compile_forms(forms, compile_args)
    rng = np.random.default_rng(1)
    dofs = rng.uniform(0.5, 1.5, V.ufl_element().dim)
    constants = np.array([1.5, -0.5])
    n = int(np.prod(shape, dtype=int))
    A = tabulate(module, compiled[0], Q.ufl_element().dim, dofs, constants)
    integral = tabulate(module, compiled[1], n, dofs, constants)
    (cell_volume,) = tabulate(module, compiled[2], 1)
    A_ref = tabulate(module, compiled[3], Q.ufl_element().dim, integral / cell_volume)
    assert np.abs(A_ref).max() > 1e-3
    np.testing.assert_allclose(A, A_ref, rtol=1e-13)


@pytest.mark.parametrize("scalar_type", [np.float64, np.complex128])
def test_cell_average_sub_kernels(scalar_type):
    """Each cell average is built by two sub-kernels, and the averages share the volume.

    The sub-kernels are cell kernels of linear forms, generated as static
    functions ahead of the kernels that call them. The proxy coefficients that
    stand for the integrals are internal to the kernels, not coefficients of
    the forms.
    """
    _, V, Q = spaces()
    u = ufl.Coefficient(V)
    v = ufl.TestFunction(Q)
    forms = [
        ufl.inner(ufl.cell_avg(u * u), v) * ufl.dx
        + ufl.inner(ufl.cell_avg(ufl.grad(u)), ufl.grad(v)) * ufl.dx,
        ufl.inner(ufl.cell_avg(u * u), v) * ufl.ds,
        ufl.inner(ufl.cell_avg(u * u), v) * ufl.dP,
    ]
    options = ffcx.options.get_options({"scalar_type": scalar_type})
    analysis = ffcx.analysis.analyze_ufl_objects(forms, options["scalar_type"])

    # Two integrals of expressions and the volume for the first form, which
    # its two integrals share, and one integral and the volume for each of the
    # others. The forms get proxies of their own, whose kernels are shared by
    # name.
    assert len(analysis.sub_kernels) == 7
    proxies = [sub_kernel.target for sub_kernel in analysis.sub_kernels]
    assert len(set(proxies)) == 7
    for sub_kernel in analysis.sub_kernels:
        assert not sub_kernel.target.is_interpolation
        assert sub_kernel.expression is None
        assert sub_kernel.form_data is not None
        assert sub_kernel.form_data.rank == 1
        (itg_data,) = sub_kernel.form_data.integral_data
        assert itg_data.integral_type == "cell"
        assert sub_kernel.target.ufl_element().embedded_superdegree == 0
    volumes = [p for p in proxies if not ffcx.analysis.required_coefficients(p.form)]
    assert len(volumes) == 3
    assert all(volume.ufl_shape == () for volume in volumes)
    numerators = [p for p in proxies if p not in volumes]
    assert sorted(p.ufl_shape for p in numerators) == [(), (), (), (2,)]
    for proxy in numerators:
        assert ffcx.analysis.sub_kernel_coefficients(proxy) == [u]

    ir = ffcx.ir.representation.compute_ir(analysis, {}, "test", options, None)
    # One kernel each for the integrals of the two expressions and for the volume
    assert len(ir.sub_kernels) == 3
    assert all(
        isinstance(sub_kernel, ffcx.ir.representation.IntegralIR) for sub_kernel in ir.sub_kernels
    )
    assert all(sub_kernel.is_sub_kernel for sub_kernel in ir.sub_kernels)
    assert sorted(sub_kernel.expression.tensor_shape for sub_kernel in ir.sub_kernels) == [
        [1],
        [1],
        [2],
    ]
    assert not any(integral.is_sub_kernel for integral in ir.integrals)
    # The cell integral reads both averages and the volume, the facet and
    # vertex integrals one and the volume
    assert sorted(len(integral.sub_kernels) for integral in ir.integrals) == [2, 2, 3]
    for integral in ir.integrals:
        assert len(integral.expression.coefficient_numbering) == 1
        assert not any(
            isinstance(c, ffcx.analysis.ProxyCoefficient)
            for c in integral.expression.coefficient_numbering
        )
        assert len({sub_kernel.target for sub_kernel in integral.sub_kernels} & set(volumes)) == 1
    for form_ir in ir.forms:
        assert form_ir.num_coefficients == 1

    (_, code), _ = ffcx.compiler.compile_ufl_objects(forms, options, namespace="test")
    sub_kernel, kernel = (
        "\nstatic void tabulate_tensor_integral_",
        "\nvoid tabulate_tensor_integral_",
    )
    assert code.count(sub_kernel) == 3
    assert code.count(kernel) == 3
    assert code.rfind(sub_kernel) < code.find(kernel)


def test_cell_average_on_facets(compile_args):
    """A cell average in a facet integral, which calls the cell sub-kernels.

    On an interior facet the average is of the cell on the side the average
    is restricted to, whose coordinates and coefficients the sub-kernels are
    passed. An average of both sides in one integrand reads each side's cell,
    and an expression restricted to both sides cannot be averaged.
    """
    mesh, V, Q = spaces()
    u = ufl.Coefficient(V)
    v = ufl.TestFunction(Q)
    f = u * u
    average = ufl.cell_avg(f)
    forms = [
        average * v * ufl.ds,
        v * ufl.ds,
        average("+") * v("+") * ufl.dS,
        v("+") * ufl.dS,
        average("-") * v("-") * ufl.dS,
        v("-") * ufl.dS,
        (average("+") - 2 * average("-")) * (v("+") + v("-")) * ufl.dS,
        (v("+") + v("-")) * ufl.dS,
        *integral_forms(f, mesh),
    ]
    compiled, module = compile_forms(forms, compile_args)
    rng = np.random.default_rng(2)
    dofs = rng.uniform(0.5, 1.5, V.ufl_element().dim)
    coords = GEOMETRIES["triangle"][2]
    entity = np.array([1, 0], dtype=np.int32)
    (integral,) = tabulate(module, compiled[8], 1, dofs, coords=coords)
    (cell_volume,) = tabulate(module, compiled[9], 1, coords=coords)

    facet = dict(coords=coords, integral_type="exterior_facet", entity=entity)
    A = tabulate(module, compiled[0], 3, dofs, **facet)
    A_ref = tabulate(module, compiled[1], 3, **facet)
    assert np.abs(A_ref).max() > 1e-3
    np.testing.assert_allclose(A, integral / cell_volume * A_ref, rtol=1e-13)

    # Two triangles of different shape sharing the facet (2, 0.3) -> (0.4, 1.7),
    # local facet 0 of both, with the same vertex order so neither side is
    # permuted. Each side has coefficients of its own, so that reading the
    # wrong cell shows.
    cells = [
        np.pad(np.array([[0.1, 0.2], [2.0, 0.3], [0.4, 1.7]]), ((0, 0), (0, 1))),
        np.pad(np.array([[2.1, 1.9], [2.0, 0.3], [0.4, 1.7]]), ((0, 0), (0, 1))),
    ]
    side_dofs = [rng.uniform(0.5, 1.5, V.ufl_element().dim) for _ in cells]
    averages = [
        tabulate(module, compiled[8], 1, w, coords=cell)[0]
        / tabulate(module, compiled[9], 1, coords=cell)[0]
        for cell, w in zip(cells, side_dofs)
    ]
    facet = dict(
        coords=np.vstack(cells),
        integral_type="interior_facet",
        entity=np.zeros(2, dtype=np.int32),
        permutation=np.zeros(2, dtype=np.uint8),
    )
    w = np.concatenate(side_dofs)
    for i, expected in [(2, averages[0]), (4, averages[1]), (6, averages[0] - 2 * averages[1])]:
        A = tabulate(module, compiled[i], 2 * 3, w, **facet)
        A_ref = tabulate(module, compiled[i + 1], 2 * 3, **facet)
        assert np.abs(A_ref).max() > 1e-3
        np.testing.assert_allclose(A, expected * A_ref, rtol=1e-13)

    with pytest.raises(ValueError, match="both sides"):
        ffcx.analysis.analyze_ufl_objects(
            [ufl.cell_avg(u("+") * u("-")) * v("+") * ufl.dS], np.float64
        )


def test_cell_average_in_expression(compile_args):
    """A cell average in an expression, which is constant over the cell."""
    mesh, V, _ = spaces()
    u = ufl.Coefficient(V)
    x = ufl.SpatialCoordinate(mesh)
    f = u * u * x[0]
    points = np.array([[0.1, 0.2], [0.5, 0.25], [0.3, 0.6]])
    compiled, module, _ = ffcx.codegeneration.jit.compile_expressions(
        [(ufl.cell_avg(f), points)], cffi_extra_compile_args=compile_args
    )
    (expression,) = compiled
    dofs = np.random.default_rng(3).uniform(0.5, 1.5, V.ufl_element().dim)
    coords = GEOMETRIES["triangle"][2]
    ffi = module.ffi
    A = np.zeros(len(points))
    expression.tabulate_tensor_float64(
        ffi.cast("double *", A.ctypes.data),
        ffi.cast("double *", dofs.ctypes.data),
        ffi.NULL,
        ffi.cast("double *", coords.flatten().ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    forms, module = compile_forms(integral_forms(f, mesh), compile_args)
    (integral,) = tabulate(module, forms[0], 1, dofs, coords=coords)
    (cell_volume,) = tabulate(module, forms[1], 1, coords=coords)
    np.testing.assert_allclose(A, integral / cell_volume, rtol=1e-13)


def test_cell_average_nested(compile_args):
    """Cell averages and interpolations nested in one another.

    An interpolation inside an averaged expression has its expression kernel
    called by the integral sub-kernel; an average inside an interpolated
    expression has its integral sub-kernels called by the expression kernel.
    Each interpolation is into the space of its field, so the identity.
    """
    _, V, Q = spaces()
    u = ufl.Coefficient(V)
    v = ufl.TestFunction(Q)
    forms = [
        ufl.cell_avg(ufl.Interpolate(u, V) * u) * v * ufl.dx,
        ufl.Interpolate(ufl.cell_avg(u * u) * u, V) * v * ufl.dx,
        ufl.cell_avg(u * u) * v * ufl.dx,
        ufl.cell_avg(u * u) * u * v * ufl.dx,
    ]
    options = ffcx.options.get_options()
    analysis = ffcx.analysis.analyze_ufl_objects(forms, options["scalar_type"])
    kinds = [sub_kernel.target.is_interpolation for sub_kernel in analysis.sub_kernels]
    # Each sub-kernel comes after those it calls: the interpolation of the
    # first form before the integral that holds it, then that form's volume,
    # and the two integrals of the second form before the interpolation that
    # holds them
    assert kinds[:3] == [True, False, False] and kinds[3:6] == [False, False, True]

    compiled, module = compile_forms(forms, compile_args)
    dofs = np.random.default_rng(4).uniform(0.5, 1.5, V.ufl_element().dim)
    A = [tabulate(module, form, Q.ufl_element().dim, dofs) for form in compiled]
    assert np.abs(A[2]).max() > 1e-3
    np.testing.assert_allclose(A[0], A[2], rtol=1e-13)
    np.testing.assert_allclose(A[1], A[3], rtol=1e-13)


def test_cell_average_of_argument_is_not_supported():
    """An expression holding an argument cannot be averaged by a sub-kernel."""
    _, V, Q = spaces()
    u, v = ufl.TrialFunction(V), ufl.TestFunction(Q)
    with pytest.raises(NotImplementedError, match="holding an argument"):
        ffcx.analysis.analyze_ufl_objects([ufl.cell_avg(u) * v * ufl.dx], np.float64)
