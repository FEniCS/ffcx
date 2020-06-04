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


@pytest.mark.parametrize("degree", [1, 2])
def test_cmap_triangle(degree, compile_args):
    """Test computation of reference coordinates for triangle cell."""
    cell = ufl.triangle
    element = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(element)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args, cache_dir=".")

    assert compiled_cmap[0].is_affine == (1 if (degree == 1) else 0)
    assert compiled_cmap[0].geometric_dimension == 2
    assert compiled_cmap[0].topological_dimension == 2

    # Reference coordinates X to basis
    phi = np.zeros(((degree + 2) * (degree + 1)) // 2, dtype=np.float64)
    phi_ptr = module.ffi.cast("double *", module.ffi.from_buffer(phi))
    X = np.array([[1 / 3, 1 / 3]], dtype=np.float64)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_cmap[0].evaluate_basis_derivatives(phi_ptr, 0, X.shape[0], X_ptr)
    assert np.isclose(sum(phi), 1.0)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[2] == 0
    assert num_entity_dofs[3] == 0

    if degree == 1:
        assert num_entity_dofs[1] == 0
    elif degree == 2:
        assert num_entity_dofs[1] == 1


@pytest.mark.parametrize("degree", [1, 2])
def test_cmap_quads(degree, compile_args):
    """Test computation of physical and reference coordinates for quadrilateral cell"""
    # Assuming FIAT Tensor Product layout of cell.

    cell = ufl.quadrilateral
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    assert compiled_cmap[0].is_affine == 0
    assert compiled_cmap[0].geometric_dimension == 2
    assert compiled_cmap[0].topological_dimension == 2

    # Reference coordinates X to basis
    phi = np.zeros((degree + 1) ** 2, dtype=np.float64)
    phi_ptr = module.ffi.cast("double *", module.ffi.from_buffer(phi))
    X = np.array([[0.5, 0.5]], dtype=np.float64)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_cmap[0].evaluate_basis_derivatives(phi_ptr, 0, X.shape[0], X_ptr)
    assert np.isclose(sum(phi), 1.0)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[3] == 0

    if degree == 1:
        assert num_entity_dofs[1] == 0
        assert num_entity_dofs[2] == 0
    elif degree == 2:
        assert num_entity_dofs[1] == 1
        assert num_entity_dofs[2] == 1


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_cmap_hex(degree, compile_args):
    """Test computation of physical and reference coordinates for hexahedron cell"""
    # Assuming FIAT Tensor Product layout of cell.

    cell = ufl.hexahedron
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    assert compiled_cmap[0].is_affine == 0
    assert compiled_cmap[0].geometric_dimension == 3
    assert compiled_cmap[0].topological_dimension == 3

    # Reference coordinates X to basis
    phi = np.zeros((degree + 1) ** 3, dtype=np.float64)
    phi_ptr = module.ffi.cast("double *", module.ffi.from_buffer(phi))
    X = np.array([[0.5, 0.5, 0.5]], dtype=np.float64)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_cmap[0].evaluate_basis_derivatives(phi_ptr, 0, X.shape[0], X_ptr)
    assert np.isclose(sum(phi), 1.0)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1

    if degree == 1:
        assert num_entity_dofs[1] == 0
        assert num_entity_dofs[2] == 0
        assert num_entity_dofs[3] == 0
    elif degree == 2:
        assert num_entity_dofs[1] == 1
        assert num_entity_dofs[2] == 1
        assert num_entity_dofs[3] == 1


@pytest.mark.parametrize("degree", [1, 2])
def test_cmap_tet(degree, compile_args):
    """Coordinate map test for tetrahedron cell"""

    cell = ufl.tetrahedron
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps(
        [mesh], cffi_extra_compile_args=compile_args)

    assert compiled_cmap[0].is_affine == (1 if (degree == 1) else 0)
    assert compiled_cmap[0].geometric_dimension == 3
    assert compiled_cmap[0].topological_dimension == 3

    # Reference coordinates X to basis
    phi = np.zeros(6 * degree - 2, dtype=np.float64)
    phi_ptr = module.ffi.cast("double *", module.ffi.from_buffer(phi))
    X = np.array([[0.25, 0.25, 0.25]], dtype=np.float64)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_cmap[0].evaluate_basis_derivatives(phi_ptr, 0, X.shape[0], X_ptr)
    assert np.isclose(sum(phi), 1.0)

    num_entity_dofs = compiled_cmap[0].create_scalar_dofmap().num_entity_dofs

    assert num_entity_dofs[0] == 1
    assert num_entity_dofs[2] == 0
    assert num_entity_dofs[3] == 0

    if degree == 1:
        assert num_entity_dofs[1] == 0
    elif degree == 2:
        assert num_entity_dofs[1] == 1
