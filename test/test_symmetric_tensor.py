# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tests for the exploit_symmetry option (ffcx.ir.integral.detect_bilinear_symmetry).

Exercised at the compiled-form level rather than by unit-testing the
detector in isolation: building a realistic BlockDataT/factorization graph
by hand would be more brittle than just compiling real forms and checking
(a) generated source for the "sym_i"/"sym_j" mirror-loop symbols the
triangular-loop path emits, and (b) that the numeric result is unaffected
except at the ULP level.
"""

import numpy as np
import pytest
import ufl
from basix.ufl import element as basix_element

import ffcx.codegeneration.jit
import ffcx.compiler
import ffcx.options

_REFERENCE_VERTICES = {
    "triangle": [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
    "tetrahedron": [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
}


def _reference_coords(cellname: str) -> np.ndarray:
    """Flattened vertex coordinates of the reference cell.

    coordinate_dofs is always stored 3D (z=0 padding for a 2D cell),
    regardless of the mesh's actual geometric dimension.
    """
    vertices = np.array(_REFERENCE_VERTICES[cellname], dtype=np.float64)
    padded = np.zeros((vertices.shape[0], 3), dtype=np.float64)
    padded[:, : vertices.shape[1]] = vertices
    return padded.flatten()


def _generated_source_contains_mirror_loop(form, options) -> bool:
    """Check whether compiling `form` with `options` emits the triangular-loop mirror."""
    opts = ffcx.options.get_options(options)
    source, _ = ffcx.compiler.compile_ufl_objects([form], options=opts)
    return "sym_i" in source[1]


def _tabulate(form, coords, exploit_symmetry, compile_args, w=None, c=None) -> np.ndarray:
    """Compile and evaluate `form`'s cell kernel, returning the reshaped local tensor."""
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [form],
        options={"scalar_type": "float64", "exploit_symmetry": exploit_symmetry},
        cffi_extra_compile_args=compile_args,
    )
    compiled_form = compiled_forms[0]
    ffi = module.ffi
    cell = module.lib.cell
    offsets = compiled_form.form_integral_offsets
    integral = compiled_form.form_integrals[offsets[cell]]

    sizes = [arg.ufl_function_space().ufl_element().dim for arg in form.arguments()]
    A = np.zeros(int(np.prod(sizes)), dtype=np.float64)
    w = np.zeros(0, dtype=np.float64) if w is None else np.ascontiguousarray(w, dtype=np.float64)
    c = np.zeros(0, dtype=np.float64) if c is None else np.ascontiguousarray(c, dtype=np.float64)
    coords = np.ascontiguousarray(coords, dtype=np.float64)

    kernel = integral.tabulate_tensor_float64
    kernel(
        ffi.cast("double *", A.ctypes.data),
        ffi.cast("double *", w.ctypes.data),
        ffi.cast("double *", c.ctypes.data),
        ffi.cast("double *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    return A.reshape(sizes)


def _poisson_form(cellname, gdim, degree):
    element = basix_element("Lagrange", cellname, degree)
    domain = ufl.Mesh(basix_element("Lagrange", cellname, 1, shape=(gdim,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    return ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx


def _elasticity_form(cellname, gdim, degree):
    element = basix_element("Lagrange", cellname, degree, shape=(gdim,))
    domain = ufl.Mesh(basix_element("Lagrange", cellname, 1, shape=(gdim,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    return ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v))) * ufl.dx


@pytest.mark.parametrize(
    "cellname,gdim,degree,form_builder",
    [
        ("triangle", 2, 1, _poisson_form),
        ("triangle", 2, 2, _poisson_form),
        ("tetrahedron", 3, 1, _poisson_form),
        ("tetrahedron", 3, 2, _poisson_form),
        ("triangle", 2, 1, _elasticity_form),
        ("triangle", 2, 2, _elasticity_form),
        ("tetrahedron", 3, 1, _elasticity_form),
        ("tetrahedron", 3, 2, _elasticity_form),
    ],
)
def test_symmetric_form_matches_full_tensor(cellname, gdim, degree, form_builder, compile_args):
    """A known-symmetric form's optimised tensor is exactly symmetric and matches the full one."""
    form = form_builder(cellname, gdim, degree)
    coords = _reference_coords(cellname)

    assert _generated_source_contains_mirror_loop(
        form, {"scalar_type": "float64", "exploit_symmetry": True}
    )

    A_full = _tabulate(form, coords, False, compile_args)
    A_sym = _tabulate(form, coords, True, compile_args)

    assert np.array_equal(A_sym, A_sym.T)
    np.testing.assert_allclose(A_sym, A_full, rtol=1e-12, atol=1e-14)


def test_coefficient_weighted_symmetric_form(compile_args):
    """A coefficient-scaled symmetric form is still recognised and matches the full tensor."""
    cellname, gdim = "triangle", 2
    element = basix_element("Lagrange", cellname, 2)
    domain = ufl.Mesh(basix_element("Lagrange", cellname, 1, shape=(gdim,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    kappa = ufl.Constant(domain, shape=(gdim, gdim))
    form = ufl.tr(kappa) * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    assert _generated_source_contains_mirror_loop(
        form, {"scalar_type": "float64", "exploit_symmetry": True}
    )

    coords = _reference_coords(cellname)
    c = np.array([2.0, 0.5, 0.5, 1.5], dtype=np.float64)  # kappa, row-major
    A_full = _tabulate(form, coords, False, compile_args, c=c)
    A_sym = _tabulate(form, coords, True, compile_args, c=c)

    assert np.array_equal(A_sym, A_sym.T)
    np.testing.assert_allclose(A_sym, A_full, rtol=1e-12, atol=1e-14)


def test_default_option_leaves_symmetry_unexploited():
    """exploit_symmetry defaults to False: no change to generated code without opting in."""
    form = _poisson_form("triangle", 2, 2)
    opts = ffcx.options.get_options({"scalar_type": "float64"})
    assert opts["exploit_symmetry"] is False
    assert not _generated_source_contains_mirror_loop(form, opts)


def test_nonsymmetric_advection_form_is_not_exploited(compile_args):
    """A genuinely non-symmetric form is correctly rejected by the detector."""
    cellname, gdim = "triangle", 2
    element = basix_element("Lagrange", cellname, 2)
    domain = ufl.Mesh(basix_element("Lagrange", cellname, 1, shape=(gdim,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    w = ufl.Constant(domain, shape=(gdim,))
    form = ufl.inner(ufl.dot(w, ufl.grad(u)), v) * ufl.dx

    assert not _generated_source_contains_mirror_loop(
        form, {"scalar_type": "float64", "exploit_symmetry": True}
    )

    coords = _reference_coords(cellname)
    c = np.array([1.0, -2.0], dtype=np.float64)
    A_full = _tabulate(form, coords, False, compile_args, c=c)
    A_sym = _tabulate(form, coords, True, compile_args, c=c)

    assert not np.allclose(A_full, A_full.T)  # sanity: the form really is asymmetric
    np.testing.assert_array_equal(A_sym, A_full)  # untouched: bit-identical, not just close


def test_mixed_test_trial_elements_not_exploited():
    """Differing test/trial elements are safely excluded (no matching swapped block exists)."""
    cellname, gdim = "triangle", 2
    domain = ufl.Mesh(basix_element("Lagrange", cellname, 1, shape=(gdim,)))
    space1 = ufl.FunctionSpace(domain, basix_element("Lagrange", cellname, 1))
    space2 = ufl.FunctionSpace(domain, basix_element("Lagrange", cellname, 2))
    u = ufl.TrialFunction(space1)
    v = ufl.TestFunction(space2)
    form = u * v * ufl.dx

    assert not _generated_source_contains_mirror_loop(
        form, {"scalar_type": "float64", "exploit_symmetry": True}
    )


def test_complex_scalar_type_not_exploited():
    """A complex-valued form is excluded (inner() conjugates one side)."""
    cellname, gdim = "triangle", 2
    element = basix_element("Lagrange", cellname, 2)
    domain = ufl.Mesh(basix_element("Lagrange", cellname, 1, shape=(gdim,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    form = ufl.inner(u, v) * ufl.dx

    assert not _generated_source_contains_mirror_loop(
        form, {"scalar_type": "complex128", "exploit_symmetry": True}
    )


def test_diagonal_part_not_exploited():
    """part='diagonal' is excluded: symmetry exploitation is only defined for the full tensor."""
    form = _poisson_form("triangle", 2, 2)
    assert not _generated_source_contains_mirror_loop(
        form, {"scalar_type": "float64", "exploit_symmetry": True, "part": "diagonal"}
    )
