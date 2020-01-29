# Copyright (C) 2018-2019 Chris Richardson and Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest

import ffcx
import ffcx.codegeneration.jit
import ufl


cells = [ufl.triangle, ufl.quadrilateral, ufl.tetrahedron, ufl.hexahedron]
degrees = [1, 2, 3, 4]


@pytest.fixture(scope="module")
def lagrange_elements():
    elements = {}
    for cell in cells:
        for degree in degrees:
            ufl_element = ufl.FiniteElement("Lagrange", cell, degree)
            compiled_elements, module = ffcx.codegeneration.jit.compile_elements([ufl_element])
            elements[(cell, degree)] = (ufl_element, compiled_elements[0], module)

    return elements


@pytest.fixture(scope="module")
def reference_points():
    points = {}

    points[ufl.triangle] = np.array([[0.0, 0.0], [0.5, 0.5], [0.25, 0.25],
                                     [1 / 3, 1 / 3], [1.0, 0.0], [0.0, 1.0]])
    points[ufl.tetrahedron] = np.array([[0.0, 0.0, 0.0], [0.5, 0.2, 0.0], [0.0, 0.0, 1.0]])
    points[ufl.quadrilateral] = np.array([[0.0, 0.0], [0.5, 0.5], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    points[ufl.hexahedron] = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [1.0, 0.0, 0.0],
                                       [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]])

    return points


@pytest.mark.parametrize("cell", cells)
@pytest.mark.parametrize("degree", degrees)
def test_dim_degree(cell, degree, lagrange_elements):
    ufl_element, compiled_element, module = lagrange_elements[(cell, degree)]

    assert compiled_element[0].geometric_dimension == cell.geometric_dimension()
    assert compiled_element[0].topological_dimension == cell.topological_dimension()
    assert ufl_element.degree() == compiled_element[0].degree


@pytest.mark.parametrize("cell", cells)
@pytest.mark.parametrize("degree", degrees)
def test_tabulate_reference_dof_coordinates(cell, degree, lagrange_elements):
    ufl_element, compiled_element, module = lagrange_elements[(cell, degree)]
    fiat_element = ffcx.fiatinterface._create_fiat_element(ufl_element)

    tdim = compiled_element[0].topological_dimension
    space_dim = compiled_element[0].space_dimension
    X = np.zeros([space_dim, tdim])
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_element[0].tabulate_reference_dof_coordinates(X_ptr)

    fiat_coordinates = np.asarray(list(sorted(L.pt_dict.keys())[0] for L in fiat_element.dual_basis()))
    assert (np.isclose(X, fiat_coordinates)).all()


@pytest.mark.parametrize("cell", cells)
@pytest.mark.parametrize("degree", degrees)
def test_evaluate_reference_basis(cell, degree, lagrange_elements, reference_points):
    ufl_element, compiled_element, module = lagrange_elements[(cell, degree)]
    fiat_element = ffcx.fiatinterface._create_fiat_element(ufl_element)

    space_dim = compiled_element[0].space_dimension
    tdim = compiled_element[0].topological_dimension

    X = reference_points[cell]
    npoint = X.shape[0]
    X_ptr = module.ffi.cast("const double *", module.ffi.from_buffer(X))
    vals = np.zeros([npoint, space_dim])
    vals_ptr = module.ffi.cast("double *", module.ffi.from_buffer(vals))

    compiled_element[0].evaluate_reference_basis(vals_ptr, npoint, X_ptr)

    fiat_vals = fiat_element.tabulate(0, X)
    assert (np.isclose(vals.T, fiat_vals[(0,) * tdim])).all()


def test_cmap_triangle():
    cell = ufl.triangle
    element = ufl.VectorElement("Lagrange", cell, 1)
    mesh = ufl.Mesh(element)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps([mesh])
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
def test_cmap_quads(degree, coords):
    # Test for first and second order quadrilateral meshes,
    # assuming FIAT Tensor Product layout of cell.

    cell = ufl.quadrilateral
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps([mesh])

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
def test_cmap_hex(degree, coords):
    # Test for first and second order quadrilateral meshes,
    # assuming FIAT Tensor Product layout of cell.

    cell = ufl.hexahedron
    e = ufl.VectorElement("Lagrange", cell, degree)
    mesh = ufl.Mesh(e)
    compiled_cmap, module = ffcx.codegeneration.jit.compile_coordinate_maps([mesh])

    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))

    x = np.array([[1 / 3, 3 / 2, 1]], dtype=np.float64)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    X = np.zeros_like(x)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_cmap[0].compute_physical_coordinates(X_ptr, X.shape[0], x_ptr, coords_ptr)

    assert(np.isclose(X[0, 0], x[0, 0]))
    assert(np.isclose(X[0, 1], 2 * x[0, 1]))
    assert(np.isclose(X[0, 2], 3 * x[0, 2]))
