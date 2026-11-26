# Copyright (C) 2026 Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
import basix
import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype

CELL = "triangle"
#: Vertices of an affine triangle that is neither axis aligned nor of unit size
COORDS = np.array([[0.1, 0.2, 0.0], [2.0, 0.3, 0.0], [0.4, 1.7, 0.0]])
#: Points in the reference triangle to evaluate expressions at
POINTS = np.array([[0.1, 0.1], [0.6, 0.2], [0.25, 0.5], [1 / 3, 1 / 3], [0.05, 0.8]])


def domain():
    return ufl.Mesh(basix.ufl.element("Lagrange", CELL, 1, shape=(2,)))


def space(mesh, family, degree, **kwargs):
    return ufl.FunctionSpace(mesh, basix.ufl.element(family, CELL, degree, **kwargs))


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


def test_unsupported_interpolations():
    """Interpolations that cannot be compiled are reported as such."""
    mesh = domain()
    V = space(mesh, "Lagrange", 1)
    u, v = ufl.Coefficient(V), ufl.TestFunction(V)

    # An expression that is only linear in the argument would need its table
    # built for each cell
    with pytest.raises(NotImplementedError, match="the expression is the argument itself"):
        ffcx.codegeneration.jit.compile_expressions([(ufl.Interpolate(u * v, V), POINTS)])

    # Including one whose argument is shaped differently from the space
    W = space(mesh, "N1curl", 1)
    with pytest.raises(NotImplementedError, match="the expression is the argument itself"):
        ffcx.codegeneration.jit.compile_expressions([(ufl.Interpolate(ufl.grad(v), W), POINTS)])

    # An expression that holds no argument would need its degrees of freedom
    # built for each cell
    with pytest.raises(NotImplementedError, match="holds no argument"):
        ffcx.codegeneration.jit.compile_expressions([(ufl.Interpolate(u * u, V), POINTS)])
