# Copyright (C) 2018-2020 Chris Richardson, Michal Habera and Jørgen S. Dokken
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import ffcx
import ffcx.codegeneration.jit
import numpy as np
import pytest

import ufl


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 4.0]], dtype=np.float64)),
                                           (2, np.array([[0, 0], [1, 0], [0, 1], [0.65, 0.65],
                                                         [-0.1, 0.5], [0.5, -0.2]], dtype=np.float64))])
def test_cmap_triangle(degree, coords, compile_args):
    """Test computation of reference coordinates for triangle cell."""
    cell = ufl.triangle
    element = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(element)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    # Reference coordinates X
    x = np.array([[1 / 3, 2 / 3]], dtype=np.float64)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    X = np.zeros_like(x)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))
    compiled_cmap[0].compute_reference_coordinates(X_ptr, X.shape[0], x_ptr, coords_ptr)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[2] == 0
    assert num_entity_dofs[3] == 0
    assert num_entity_dofs[1] == degree - 1

    if degree == 1:
        assert np.isclose(X[0, 0], 1 / 6)
        assert np.isclose(X[0, 1], 1 / 6)

    # Convert back to reference coordinates
    Y = np.zeros_like(X)
    Y_ptr = module.ffi.cast("double *", module.ffi.from_buffer(Y))
    retcode = compiled_cmap[0].compute_reference_coordinates(Y_ptr, Y.shape[0], x_ptr, coords_ptr)
    assert np.isclose(X, Y).all()
    assert retcode == 0


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0, 0], [3, 0], [0, 2], [3, 2]], dtype=np.float64)),
                                           (2, np.array([[0, 0], [3, 0], [1.5, 0], [0, 2],
                                                         [3, 2], [1.5, 2], [0, 1], [3, 1], [1.5, 1]],
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


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0, 0], [3, 0], [0, 2], [3.1, 2.1]], dtype=np.float64)),
                                           (2, np.array([[0, 0], [3, 0], [1.5, 0], [0, 2],
                                                         [3.1, 2.1], [1.5, 2], [0, 1], [3, 1], [1.5, 1]],
                                                        dtype=np.float64))])
def test_cmap_quad_distorted(degree, coords, compile_args):
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


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0, 0, 0], [1, 0, 0],
                                                         [0, 2, 0], [1, 2, 0],
                                                         [0, 0, 3], [1, 0, 3],
                                                         [0, 2, 3], [1.1, 2.2, 3.3]], dtype=np.float64)),
                                           (2, np.array([[0, 0, 0], [0, 0, 1.1], [-0.1, 0, 0.5],
                                                         [0, 0.9, 0.1], [0, 0.8, 1], [-0.2, 1, 0.5],
                                                         [0, 0.5, -0.1], [-0.1, 0.3, 1.2], [-0.3, 0.5, 0.5],
                                                         [1, 0, 0], [1, 0, 1.1], [1, 0, 0.5],
                                                         [1, 0.9, 0.1], [1, 0.8, 1], [1, 1, 0.5],
                                                         [1, 0.5, -0.1], [1, 0.3, 1.4], [1, 0.5, 0.5],
                                                         [0.5, 0, 0], [0.5, 0, 1.1], [0.5, 0, 0.5],
                                                         [0.5, 0.9, 0.1], [0.5, 0.8, 1], [0.5, 1, 0.5],
                                                         [0.5, 0.5, -0.1], [0.5, 0.3, 1.4], [0.5, 0.5, 0.5],
                                                         ], dtype=np.float64))])
def test_cmap_hex_distorted(degree, coords, compile_args):
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

    # Convert back to reference coordinates
    Y = np.zeros_like(X)
    Y_ptr = module.ffi.cast("double *", module.ffi.from_buffer(Y))
    retcode = compiled_cmap[0].compute_reference_coordinates(Y_ptr, Y.shape[0], x_ptr, coords_ptr)
    assert np.isclose(X, Y).all()
    assert retcode == 0


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0, 0, 0], [1, 0, 0],
                                                         [0, 2, 0], [0, 0, 3]],
                                                        dtype=np.float64)),
                                           (2, np.array([[0, 0, 0], [1, 0, 0],
                                                         [0, 2, 0], [0, 0, 3],
                                                         [0, 1, 1.5], [0.5, 0, 1.5],
                                                         [0.5, 1, 0], [0, 0, 1.5],
                                                         [0, 1, 0], [0.5, 0, 0]],
                                                        dtype=np.float64))])
def test_cmap_tet(degree, coords, compile_args):
    """Test computation of physical and reference coordinates for hexahedron cell"""
    # Assuming FIAT Tensor Product layout of cell.

    cell = ufl.tetrahedron
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))

    # Reference coordinates X
    X = np.array([[1 / 3, 1 / 2, 1]], dtype=np.float64)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    # Physical coordinates x
    x = np.zeros_like(X)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    compiled_cmap[0].compute_physical_coordinates(x_ptr, x.shape[0], X_ptr, coords_ptr)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[1] == max(0, degree - 1)
    assert num_entity_dofs[2] == max(0, degree - 2)
    assert num_entity_dofs[3] == max(0, (degree - 2) * (degree - 1) / 2)

    assert(np.isclose(x[0, 0], X[0, 0]))
    assert(np.isclose(x[0, 1], 2 * X[0, 1]))
    assert(np.isclose(x[0, 2], 3 * X[0, 2]))

    # Convert back to reference coordinates
    Y = np.zeros_like(X)
    Y_ptr = module.ffi.cast("double *", module.ffi.from_buffer(Y))
    retcode = compiled_cmap[0].compute_reference_coordinates(Y_ptr, Y.shape[0], x_ptr, coords_ptr)
    assert np.isclose(X, Y).all()
    assert retcode == 0


@pytest.mark.parametrize("degree,coords", [(1, np.array([[0.1, -0.2, 0], [1.1, 0, 0],
                                                         [0, 1.7, 0.1], [0.1, 0.2, 2]],
                                                        dtype=np.float64)),
                                           (2, np.array([[0.0, 0.0, 0.0], [1, 0, 0],
                                                         [0, 1, 0], [0, 0, 1],
                                                         [0, 0.55, 0.55], [0.65, -0.1, 0.53],
                                                         [0.51, 0.54, -0.2], [-0.1, -0.05, 0.52],
                                                         [0, 0.5, 0], [0.5, -0.5, -0.1]]))])
def test_cmap_tet_distorted(degree, coords, compile_args):
    """Test computation of physical and reference coordinates for tetrahedron cell"""
    # Assuming FIAT Tensor Product layout of cell.

    cell = ufl.tetrahedron
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))

    # Reference coordinates X
    X = np.array([[1 / 3, 1 / 2, 1]], dtype=np.float64)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    # Physical coordinates x
    x = np.zeros_like(X)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    compiled_cmap[0].compute_physical_coordinates(x_ptr, x.shape[0], X_ptr, coords_ptr)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[1] == max(0, degree - 1)
    assert num_entity_dofs[2] == max(0, degree - 2)
    assert num_entity_dofs[3] == max(0, (degree - 2) * (degree - 1) / 2)

    # Convert back to reference coordinates
    Y = np.zeros_like(X)
    Y_ptr = module.ffi.cast("double *", module.ffi.from_buffer(Y))
    retcode = compiled_cmap[0].compute_reference_coordinates(Y_ptr, Y.shape[0], x_ptr, coords_ptr)
    assert np.isclose(X, Y).all()
    assert retcode == 0
