# Copyright (C) 2018-2020 Chris Richardson and Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np

import ffcx
import ffcx.codegeneration.jit
import pytest
import ufl


def test_cmap_triangle(compile_args):
    """Test computation of reference coordinates for triangle cell."""
    cell = ufl.triangle
    element = ufl.VectorElement("Lagrange", cell, 1)
    mesh = ufl.Mesh(element)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    # Reference coordinates X
    x = np.array([[0.5, 0.5]], dtype=np.float64)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    X = np.zeros_like(x)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    coords = np.array([0.0, 0.0, 2.0, 0.0, 0.0, 4.0], dtype=np.float64)
    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))
    compiled_cmap[0].compute_reference_coordinates(X_ptr, X.shape[0], x_ptr, coords_ptr)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[1] == 0
    assert num_entity_dofs[2] == 0
    assert num_entity_dofs[3] == 0

    assert np.isclose(X[0, 0], 0.25)
    assert np.isclose(X[0, 1], 0.125)


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0, 0], [0, 2], [3, 0], [3, 2]], dtype=np.float64)),
                                           (2, np.array([[0, 0], [0, 2], [0, 1], [3, 0],
                                                         [3, 2], [3, 1], [1.5, 0], [1.5, 2], [1.5, 1]],
                                                        dtype=np.float64))])
def test_cmap_quads(degree, coords, compile_args):
    """Test computation of physical and reference coordinates for quadrilateral cell"""
    # Assuming FIAT Tensor Product layout of cell.

    cell = ufl.quadrilateral
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))

    # Reference coordinates X
    X = np.array([[1 / 3, 1 / 3]], dtype=np.float64)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    # Physical coordinates x
    x = np.zeros_like(X)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))

    compiled_cmap[0].compute_physical_coordinates(x_ptr, x.shape[0], X_ptr, coords_ptr)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[1] == degree - 1
    assert num_entity_dofs[2] == (degree - 1) ** 2
    assert num_entity_dofs[3] == 0

    assert(np.isclose(x[0, 0], 3 * X[0, 0]))
    assert(np.isclose(x[0, 1], 2 * X[0, 1]))

    # Convert back to reference coordinates
    Y = np.zeros_like(X)
    Y_ptr = module.ffi.cast("double *", module.ffi.from_buffer(Y))
    retcode = compiled_cmap[0].compute_reference_coordinates(Y_ptr, Y.shape[0], x_ptr, coords_ptr)
    assert np.isclose(X, Y).all()
    assert retcode == 0


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
    """Test computation of physical and reference coordinates for hexahedron cell"""
    # Assuming FIAT Tensor Product layout of cell.

    cell = ufl.hexahedron
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))

    # Reference coordinates X
    X = np.array([[1 / 3, 3 / 2, 1]], dtype=np.float64)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    # Physical coordinates x
    x = np.zeros_like(X)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    compiled_cmap[0].compute_physical_coordinates(x_ptr, x.shape[0], X_ptr, coords_ptr)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[1] == degree - 1
    assert num_entity_dofs[2] == (degree - 1) ** 2
    assert num_entity_dofs[3] == (degree - 1) ** 3

    assert(np.isclose(x[0, 0], X[0, 0]))
    assert(np.isclose(x[0, 1], 2 * X[0, 1]))
    assert(np.isclose(x[0, 2], 3 * X[0, 2]))

    # Convert back to reference coordinates
    Y = np.zeros_like(X)
    Y_ptr = module.ffi.cast("double *", module.ffi.from_buffer(Y))
    retcode = compiled_cmap[0].compute_reference_coordinates(Y_ptr, Y.shape[0], x_ptr, coords_ptr)
    assert np.isclose(X, Y).all()
    assert retcode == 0
