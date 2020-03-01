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

families = ["Lagrange", "Brezzi-Douglas-Marini", "Raviart-Thomas", "N1curl", "N2curl"]
cells = [ufl.triangle, ufl.quadrilateral, ufl.tetrahedron, ufl.hexahedron]
degrees = [1, 2, 3]


@pytest.fixture(scope="module")
def reference_points():
    """Returns an example reference points in reference cells"""
    points = {}

    points[ufl.interval] = np.array([[0.0], [0.3], [0.9]])
    points[ufl.triangle] = np.array([[0.0, 0.0], [0.5, 0.5], [0.25, 0.25],
                                     [1 / 3, 1 / 3], [1.0, 0.0], [0.0, 1.0]])
    points[ufl.tetrahedron] = np.array([[0.0, 0.0, 0.0], [0.5, 0.2, 0.0], [0.0, 0.0, 1.0]])
    points[ufl.quadrilateral] = np.array([[0.0, 0.0], [0.5, 0.5], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    points[ufl.hexahedron] = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [1.0, 0.0, 0.0],
                                       [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]])

    return points


def test_dim_degree(compiled_element):
    ufl_element, compiled_element, module = compiled_element
    cell = ufl_element.cell()

    assert compiled_element[0].geometric_dimension == cell.geometric_dimension()
    assert compiled_element[0].topological_dimension == cell.topological_dimension()
    assert ufl_element.degree() == compiled_element[0].degree


def test_tabulate_reference_dof_coordinates(compiled_element):
    ufl_element, compiled_element, module = compiled_element

    if ufl_element.family() != "Lagrange":
        pytest.skip("Cannot tabulate dofs for non-lagrange FE.")

    fiat_element = ffcx.fiatinterface._create_fiat_element(ufl_element)

    tdim = compiled_element[0].topological_dimension
    space_dim = compiled_element[0].space_dimension
    X = np.zeros([space_dim, tdim])
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_element[0].tabulate_reference_dof_coordinates(X_ptr)

    fiat_coordinates = np.asarray(list(sorted(L.pt_dict.keys())[0] for L in fiat_element.dual_basis()))
    assert (np.isclose(X, fiat_coordinates)).all()


def test_evaluate_reference_basis(compiled_element, reference_points):
    ufl_element, compiled_element, module = compiled_element

    fiat_element = ffcx.fiatinterface._create_fiat_element(ufl_element)

    space_dim = compiled_element[0].space_dimension
    tdim = compiled_element[0].topological_dimension

    # For vector/tensor valued basis this is not 1
    value_size = np.product(fiat_element.value_shape(), dtype=np.int)

    X = reference_points[ufl_element.cell()]
    npoint = X.shape[0]
    X_ptr = module.ffi.cast("const double *", module.ffi.from_buffer(X))
    vals = np.zeros([npoint, space_dim, value_size])
    vals_ptr = module.ffi.cast("double *", module.ffi.from_buffer(vals))

    compiled_element[0].evaluate_reference_basis(vals_ptr, npoint, X_ptr)

    fiat_vals = fiat_element.tabulate(0, X)

    # FFC does some reordering and slicing wrt. FIAT
    vals = np.transpose(vals, axes=[1, 2, 0])
    if value_size == 1:
        vals = vals[:, 0, :]

    assert (np.isclose(vals, fiat_vals[(0,) * tdim])).all()
