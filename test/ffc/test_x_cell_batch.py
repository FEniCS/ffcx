# -*- coding: utf-8 -*-
# Copyright (C) 2018 Fabian Löschner
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Test for cell batching or vectorization features"""

import numpy
import ctypes
import cffi
from copy import deepcopy

import ufl
from ufl import inner, dot, grad, tr, det, ln, dx, ds

import ffc.backends.ufc.jit
from ffc.fiatinterface import create_element

# FIXME: Test multiple batch sizes (currently a batch size of 4 cells is hardcoded)
# FIXME: Correctly test forms with multiple tabulate_tensor functions
#        (currently only uses "create_default_cell_integral()")


def _compute_shapes(form: ufl.Form):
    """Computes the shapes of the cell tensor, coefficient & coordinate arrays for a form"""

    # Cell tensor
    elements = tuple(arg.ufl_element() for arg in form.arguments())
    fiat_elements = (create_element(element) for element in elements)
    element_dims = tuple(fe.space_dimension() for fe in fiat_elements)
    A_shape = element_dims

    # Coefficient arrays
    ws_shapes = []
    for coefficient in form.coefficients():
        w_element = coefficient.ufl_element()
        w_dim = create_element(w_element).space_dimension()
        ws_shapes.append(w_dim)

    # Coordinate dofs array
    num_vertices_per_cell = form.ufl_cell().num_vertices()
    gdim = form.ufl_cell().geometric_dimension()
    coords_shape = (num_vertices_per_cell, gdim)

    return A_shape, ws_shapes, coords_shape


def _generate_data(n: int, A_shape, ws_shapes, coords_shape):
    """Generates n random test data sets (with random values) of the specified shapes"""

    A_cells = []
    ws_cells = []
    coords_cells = []

    for i in range(n):
        A = numpy.empty(A_shape, dtype=numpy.float64)
        A[:] = numpy.random.random_sample(A.shape)
        A_cells.append(A)

        ws = []
        for w_shape in ws_shapes:
            w = numpy.empty(w_shape, dtype=numpy.float64)
            w.fill(1.0)
            ws.append(w)
        ws_cells.append(ws)

        coords = numpy.empty(coords_shape, dtype=numpy.float64)
        coords[:] = numpy.random.random_sample(coords.shape)
        coords_cells.append(coords)

    return A_cells, ws_cells, coords_cells


def _make_interleaved(A_cells):
    """Transforms a list of arrays for different cells to a single interleaved array"""

    # Handle lists of coefficients arrays
    if isinstance(A_cells[0], list):
        A = []
        for a in zip(*A_cells):
            A.append(_make_interleaved(list(a)))
    # Handle cell tensors and coord arrays
    else:
        A = numpy.empty(A_cells[0].shape + (len(A_cells),), dtype=numpy.float64)
        for i in range(len(A_cells)):
            A[..., i] = A_cells[i]

    return A


def _test_form(compiled_form, A, ws, coords):
    """Runs a compiled form with the supplied arguments"""

    # Get pointers for tabulate_tensor arguments
    double_ptr_t = ctypes.POINTER(ctypes.c_double)

    _A_ptr = A.ctypes.data
    _coords_ptr = coords.ctypes.data

    # Create array of pointers to coefficient arrays
    _ws = (double_ptr_t * len(ws))()
    for i, w in enumerate(ws):
        _ws[i] = ctypes.cast(w.ctypes.data, double_ptr_t)

    _ws_ptr = ctypes.addressof(_ws)

    # Get tabulate_tensor
    cell_integral = compiled_form.create_default_cell_integral()
    tabulate_tensor = cell_integral.tabulate_tensor

    ffi = cffi.FFI()
    tabulate_tensor(ffi.cast("double*", _A_ptr),
                    ffi.cast("double**", _ws_ptr),
                    ffi.cast("double*", _coords_ptr),
                    ffi.cast("int", 1))


def _test_runner(forms, names):
    n = 4

    print("Compiling reference...")
    ref_compiled_forms, _ = ffc.backends.ufc.jit.compile_forms(forms)
    print("Compiling vectorized...")
    vec_compiled_forms, _ = ffc.backends.ufc.jit.compile_forms(forms, parameters={"cell_batch_size": n})
    print("Compiling vector exts...")
    ext_compiled_forms, _ = ffc.backends.ufc.jit.compile_forms(forms, parameters={"cell_batch_size": n,
                                                                                  "enable_cross_cell_gcc_ext": 1})

    for i, (form, name) in enumerate(zip(forms, names)):
        shapes = _compute_shapes(form)
        A_cells, ws_cells, coords_cells = _generate_data(n, *shapes)

        # Reference code
        reference_compiled = ref_compiled_forms[i]
        As1 = deepcopy(A_cells)
        ws1 = deepcopy(ws_cells)
        coords1 = deepcopy(coords_cells)

        for j in range(n):
            _test_form(reference_compiled, As1[j], ws1[j], coords1[j])

        A1 = _make_interleaved(As1)
        A1_norm = numpy.linalg.norm(A1)

        # for-loop batched code
        vec_compiled = vec_compiled_forms[i]
        A2 = _make_interleaved(A_cells)
        ws2 = _make_interleaved(ws_cells)
        coords2 = _make_interleaved(coords_cells)

        _test_form(vec_compiled, A2, ws2, coords2)

        A2_norm = numpy.linalg.norm(A2)
        assert numpy.isclose(A1_norm, A2_norm), "{}: norm wrong " \
                                                "(Form {}, cell_batch_size={}, no gcc-exts)".format(name, form, n)
        assert numpy.allclose(A1, A2), "{}: tensor entries wrong " \
                                       "(Form {}, cell_batch_size={}, no gcc-exts)".format(name, form, n)

        # vector extension batched code
        ext_compiled = ext_compiled_forms[i]
        A3 = _make_interleaved(A_cells)
        ws3 = _make_interleaved(ws_cells)
        coords3 = _make_interleaved(coords_cells)

        _test_form(ext_compiled, A3, ws3, coords3)

        A3_norm = numpy.linalg.norm(A3)
        assert numpy.isclose(A1_norm, A3_norm), "{}: norm wrong " \
                                                "(Form {}, cell_batch_size={}, gcc-exts)".format(name, form, n)
        assert numpy.allclose(A1, A3), "{}: tensor entries wrong " \
                                       "(Form {}, cell_batch_size={}, gcc-exts)".format(name, form, n)


def test_poisson():
    cell = ufl.tetrahedron
    element_p1 = ufl.FiniteElement("Lagrange", cell, 1)
    element_p2 = ufl.FiniteElement("Lagrange", cell, 2)

    u = ufl.TrialFunction(element_p2)
    v = ufl.TestFunction(element_p2)

    c = ufl.Coefficient(element_p1)
    f = ufl.Coefficient(element_p2)

    a = inner(c * grad(u), grad(v)) * dx
    L = f * v * dx

    _test_runner([a, L], ["Poisson bilinear form", "Poisson linear form"])


def test_hyperelasticity():
    cell = ufl.tetrahedron
    vector = ufl.VectorElement("Lagrange", cell, 1)
    scalar = ufl.FiniteElement("Lagrange", cell, 1)

    # Coefficients
    v = ufl.TestFunction(vector)  # Test function
    du = ufl.TrialFunction(vector)  # Incremental displacement
    u = ufl.Coefficient(vector)  # Displacement from previous iteration

    B = ufl.Coefficient(vector)  # Body force per unit volume
    T = ufl.Coefficient(vector)  # Traction force on the boundary

    # Kinematics
    d = u.geometric_dimension()
    F = ufl.Identity(d) + grad(u)  # Deformation gradient
    C = F.T * F  # Right Cauchy-Green tensor

    # Invariants of deformation tensors
    Ic = tr(C)
    J = det(F)

    # Elasticity parameters
    mu, lmbda = ufl.Coefficient(scalar), ufl.Coefficient(scalar)

    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J)) ** 2

    # Total potential energy
    Pi = psi * dx - dot(B, u) * dx - dot(T, u) * ds

    # Compute first variation of Pi (directional derivative about u in the direction of v)
    F = ufl.derivative(Pi, u, v)

    # Compute Jacobian of F
    J = ufl.derivative(F, u, du)

    _test_runner([J, F], ["Hyperelasticity Jacobian", "Hyperelasticity forces"])
