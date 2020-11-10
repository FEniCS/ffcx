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
