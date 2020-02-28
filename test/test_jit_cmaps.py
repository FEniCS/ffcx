# Copyright (C) 2018-2020 Chris Richardson and Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest

import ffcx
import ffcx.codegeneration.jit
import ufl


def test_cmap_triangle(compile_args):
    cell = ufl.triangle
    element = ufl.VectorElement("Lagrange", cell, 1)
    mesh = ufl.Mesh(element)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)
    x = np.array([[0.5, 0.5]], dtype=np.float64)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    X = np.zeros_like(x)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    coords = np.array([0.0, 0.0, 2.0, 0.0, 0.0, 4.0], dtype=np.float64)
    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))
    compiled_cmap[0].compute_reference_coordinates(X_ptr, X.shape[0], x_ptr, coords_ptr)

    assert(np.isclose(X[0, 0], 0.25))
    assert(np.isclose(X[0, 1], 0.125))


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0, 0], [0, 2], [3, 0], [3, 2]], dtype=np.float64)),
                                           (2, np.array([[0, 0], [0, 2], [0, 1], [3, 0],
                                                         [3, 2], [3, 1], [1.5, 0], [1.5, 2], [1.5, 1]],
                                                        dtype=np.float64))])
def test_cmap_quads(degree, coords, compile_args):
    # Test for first and second order quadrilateral meshes,
    # assuming FIAT Tensor Product layout of cell.

    cell = ufl.quadrilateral
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))

    x = np.array([[1 / 3, 1 / 3]], dtype=np.float64)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    X = np.zeros_like(x)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))

    compiled_cmap[0].compute_physical_coordinates(X_ptr, X.shape[0], x_ptr, coords_ptr)

    assert(np.isclose(X[0, 0], 3 * x[0, 0]))
    assert(np.isclose(X[0, 1], 2 * x[0, 1]))


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0, 0, 0], [0, 0, 3],
                                                         [0, 2, 0], [0, 2, 3],
                                                         [1, 0, 0], [1, 0, 3],
                                                         [1, 2, 0], [1, 2, 3]], dtype=np.float64)),
                                           (2, np.array([[0, 0, 0], [0, 0, 3], [0, 0, 1.5],
                                                         [0, 2, 0], [0, 2, 3], [0, 2, 1.5],
                                                         [0, 1, 0], [0, 1, 3], [0, 1, 1.5],
                                                         [1, 0, 0], [1, 0, 3], [1, 0, 1.5],
                                                         [1, 2, 0], [1, 2, 3], [1, 2, 1.5],
                                                         [1, 1, 0], [1, 1, 3], [1, 1, 1.5],
                                                         [0.5, 0, 0], [0.5, 0, 3], [0.5, 0, 1.5],
                                                         [0.5, 2, 0], [0.5, 2, 3], [0.5, 2, 1.5],
                                                         [0.5, 1, 0], [0.5, 1, 3], [0.5, 1, 1.5]], dtype=np.float64))])
def test_cmap_hex(degree, coords, compile_args):
    # Test for first and second order quadrilateral meshes,
    # assuming FIAT Tensor Product layout of cell.

    cell = ufl.hexahedron
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))

    x = np.array([[1 / 3, 3 / 2, 1]], dtype=np.float64)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    X = np.zeros_like(x)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_cmap[0].compute_physical_coordinates(X_ptr, X.shape[0], x_ptr, coords_ptr)

    assert(np.isclose(X[0, 0], x[0, 0]))
    assert(np.isclose(X[0, 1], 2 * x[0, 1]))
    assert(np.isclose(X[0, 2], 3 * x[0, 2]))
