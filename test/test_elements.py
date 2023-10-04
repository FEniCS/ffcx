# Copyright (C) 2007-2017 Anders Logg and Garth N. Wells
#
# This file is part of FFCx.
#
# FFCx is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFCx is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFCx. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Marie E. Rognes, 2010
# Modified by Lizao Li, 2016
"Unit tests for FFCx"


import numpy as np
import pytest

import basix.ufl


def element_coords(cell):
    if cell == "interval":
        return [(0,), (1,)]
    elif cell == "triangle":
        return [(0, 0), (1, 0), (0, 1)]
    elif cell == "tetrahedron":
        return [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    elif cell == "quadrilateral":
        return [(0, 0), (1, 0), (0, 1), (1, 1)]
    elif cell == "hexahedron":
        return [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
    else:
        RuntimeError("Unknown cell type")


def random_point(shape):
    w = np.random.random(len(shape))
    return sum([np.array(shape[i]) * w[i] for i in range(len(shape))]) / sum(w)


@pytest.mark.parametrize("degree, expected_dim", [(1, 3), (2, 6), (3, 10)])
def test_continuous_lagrange(degree, expected_dim):
    "Test space dimensions of continuous Lagrange elements."
    P = basix.ufl.element("Lagrange", "triangle", degree)
    assert P.dim == expected_dim


@pytest.mark.parametrize("degree, expected_dim", [(1, 4), (2, 9), (3, 16)])
def xtest_continuous_lagrange_quadrilateral(degree, expected_dim):
    "Test space dimensions of continuous TensorProduct elements (quadrilateral)."
    P = basix.ufl.element("Lagrange", "quadrilateral", degree)
    assert P.dim == expected_dim


@pytest.mark.parametrize("degree, expected_dim", [(1, 4), (2, 9), (3, 16)])
def xtest_continuous_lagrange_quadrilateral_spectral(degree, expected_dim):
    "Test space dimensions of continuous TensorProduct elements (quadrilateral)."
    P = basix.ufl.element("Lagrange", "quadrilateral", degree, variant="spectral")
    assert P.dim == expected_dim


@pytest.mark.parametrize("degree, expected_dim", [(0, 1), (1, 3), (2, 6), (3, 10)])
def test_discontinuous_lagrange(degree, expected_dim):
    "Test space dimensions of discontinuous Lagrange elements."
    P = basix.ufl.element("DG", "triangle", degree)
    assert P.dim == expected_dim


@pytest.mark.parametrize("degree, expected_dim",
                         [(0, 3), (1, 9), (2, 18), (3, 30)])
def test_regge(degree, expected_dim):
    "Test space dimensions of generalized Regge element."
    P = basix.ufl.element("Regge", "triangle", degree)
    assert P.dim == expected_dim


@pytest.mark.parametrize("degree, expected_dim",
                         [(0, 3), (1, 9), (2, 18), (3, 30)])
def xtest_hhj(degree, expected_dim):
    "Test space dimensions of Hellan-Herrmann-Johnson element."
    P = basix.ufl.element("HHJ", "triangle", degree)
    assert P.dim == expected_dim


class TestFunctionValues():
    """These tests examine tabulate gives the correct answers for a the
supported (non-mixed) for low degrees"""

    # FIXME: Add tests for NED and BDM/RT in 3D.

    # Shape (basis) functions on reference element
    reference_interval_1 = [lambda x: 1 - x[0], lambda x: x[0]]
    reference_triangle_1 = [lambda x: 1 - x[0] - x[1], lambda x: x[0], lambda x: x[1]]
    reference_tetrahedron_1 = [lambda x: 1 - x[0] - x[1] - x[2], lambda x: x[0],
                               lambda x: x[1], lambda x: x[2]]
    reference_triangle_bdm1 = [lambda x: (2 * x[0], -x[1]),
                               lambda x: (-x[0], 2 * x[1]),
                               lambda x: (2 - 2 * x[0] - 3 * x[1], x[1]),
                               lambda x: (- 1 + x[0] + 3 * x[1], - 2 * x[1]),
                               lambda x: (-x[0], -2 + 3 * x[0] + 2 * x[1]),
                               lambda x: (2 * x[0], 1 - 3 * x[0] - x[1])]
    reference_triangle_rt1 = [lambda x: (-x[0], -x[1]), lambda x: (x[0] - 1, x[1]),
                              lambda x: (-x[0], 1 - x[1])]
    reference_triangle_rt2 = [lambda x: (x[0] - 3 * x[0]**2, x[1] - 3 * x[0] * x[1]),
                              lambda x: (x[0] - 3 * x[0] * x[1], x[1] - 3 * x[1]**2),
                              lambda x: (-2 + 5 * x[0] + 3 * x[1] - 3 * x[0] * x[1] - 3 * x[0]**2,
                                         2 * x[1] - 3 * x[0] * x[1] - 3 * x[1]**2),
                              lambda x: (1.0 - x[0] - 3 * x[1] + 3 * x[0] * x[1], x[1] + 3 * x[1]**2),
                              lambda x: (-2 * x[0] + 3 * x[0] * x[1] + 3 * x[0] ** 2,
                                         2 - 3 * x[0] - 5 * x[1] + 3 * x[0] * x[1] + 3 * x[1]**2),
                              lambda x: (x[0] - 3 * x[0]**2,
                                         -1 + 3 * x[0] + x[1] - 3 * x[0] * x[1]),
                              lambda x: (-6 * x[0] + 3 * x[0] * x[1] + 6 * x[0]**2,
                                         -3 * x[1] + 6 * x[0] * x[1] + 3 * x[1]**2),
                              lambda x: (-3 * x[0] + 6 * x[0] * x[1] + 3 * x[0]**2,
                                         -6 * x[1] + 3 * x[0] * x[1] + 6 * x[1]**2)]
    reference_triangle_ned1 = [lambda x: (-x[1], x[0]), lambda x: (x[1], 1 - x[0]),
                               lambda x: (1.0 - x[1], x[0])]
    reference_tetrahedron_rt1 = [lambda x: (2 ** 0.5 * x[0], 2 ** 0.5 * x[1], 2 ** 0.5 * x[2]),
                                 lambda x: (2 ** 0.5 - 2 ** 0.5 * x[0], -2 ** 0.5 * x[1], -2 ** 0.5 * x[2]),
                                 lambda x: (2 ** 0.5 * x[0], 2 ** 0.5 * x[1] - 2 ** 0.5, 2 ** 0.5 * x[2]),
                                 lambda x: (-2 ** 0.5 * x[0], -2 ** 0.5 * x[1], 2 ** 0.5 - 2 ** 0.5 * x[2])]
    reference_tetrahedron_bdm1 = [lambda x: (-3 * x[0], x[1], x[2]),
                                  lambda x: (x[0], -3 * x[1], x[2]),
                                  lambda x: (x[0], x[1], -3 * x[2]),
                                  lambda x: (-3.0 + 3 * x[0] + 4 * x[1] + 4 * x[2], -x[1], -x[2]),
                                  lambda x: (1.0 - x[0] - 4 * x[1], 3 * x[1], -x[2]),
                                  lambda x: (1.0 - x[0] - 4 * x[2], -x[1], 3 * x[2]),
                                  lambda x: (x[0], 3.0 - 4 * x[0] - 3 * x[1] - 4 * x[2], x[2]),
                                  lambda x: (-3 * x[0], -1.0 + 4 * x[0] + x[1], x[2]),
                                  lambda x: (x[0], -1.0 + x[1] + 4 * x[2], -3 * x[2]),
                                  lambda x: (-x[0], -x[1], -3.0 + 4 * x[0] + 4 * x[1] + 3 * x[2]),
                                  lambda x: (3 * x[0], -x[1], 1.0 - 4 * x[0] - x[2]),
                                  lambda x: (-x[0], 3 * x[1], 1.0 - 4 * x[1] - x[2])]
    reference_tetrahedron_ned1 = [lambda x: (0.0, -x[2], x[1]),
                                  lambda x: (-x[2], 0.0, x[0]),
                                  lambda x: (-x[1], x[0], 0.0),
                                  lambda x: (x[2], x[2], 1.0 - x[0] - x[1]),
                                  lambda x: (x[1], 1.0 - x[0] - x[2], x[1]),
                                  lambda x: (1.0 - x[1] - x[2], x[0], x[0])]
    reference_quadrilateral_1 = [lambda x: (1 - x[0]) * (1 - x[1]),
                                 lambda x: (1 - x[0]) * x[1],
                                 lambda x: x[0] * (1 - x[1]),
                                 lambda x: x[0] * x[1]]
    reference_hexahedron_1 = [lambda x: (1 - x[0]) * (1 - x[1]) * (1 - x[2]),
                              lambda x: (1 - x[0]) * (1 - x[1]) * x[2],
                              lambda x: (1 - x[0]) * x[1] * (1 - x[2]),
                              lambda x: (1 - x[0]) * x[1] * x[2],
                              lambda x: x[0] * (1 - x[1]) * (1 - x[2]),
                              lambda x: x[0] * (1 - x[1]) * x[2],
                              lambda x: x[0] * x[1] * (1 - x[2]),
                              lambda x: x[0] * x[1] * x[2]]

    # Tests to perform
    tests = [("Lagrange", "interval", 1, reference_interval_1),
             ("Lagrange", "triangle", 1, reference_triangle_1),
             ("Lagrange", "tetrahedron", 1, reference_tetrahedron_1),
             # ("Lagrange", "quadrilateral", 1, reference_quadrilateral_1),
             # ("Lagrange", "hexahedron", 1, reference_hexahedron_1),
             ("Discontinuous Lagrange", "interval", 1, reference_interval_1),
             ("Discontinuous Lagrange", "triangle", 1, reference_triangle_1),
             ("Discontinuous Lagrange", "tetrahedron", 1, reference_tetrahedron_1),
             # ("Brezzi-Douglas-Marini", "triangle", 1, reference_triangle_bdm1),
             ("Raviart-Thomas", "triangle", 1, reference_triangle_rt1),
             # ("Raviart-Thomas", "triangle", 2, reference_triangle_rt2),
             # ("Discontinuous Raviart-Thomas", "triangle", 1, reference_triangle_rt1),
             # ("Discontinuous Raviart-Thomas", "triangle", 2, reference_triangle_rt2),
             ("N1curl", "triangle", 1, reference_triangle_ned1),
             ("Raviart-Thomas", "tetrahedron", 1, reference_tetrahedron_rt1),
             # ("Discontinuous Raviart-Thomas", "tetrahedron", 1, reference_tetrahedron_rt1),
             # ("Brezzi-Douglas-Marini", "tetrahedron", 1, reference_tetrahedron_bdm1),
             ("N1curl", "tetrahedron", 1, reference_tetrahedron_ned1)]

    @pytest.mark.parametrize("family, cell, degree, reference", tests)
    def test_values(self, family, cell, degree, reference):
        # Create element
        e = basix.ufl.element(family, cell, degree)

        # Get some points and check basis function values at points
        points = [random_point(element_coords(cell)) for i in range(5)]
        for x in points:
            table = e.tabulate(0, np.array([x], dtype=np.float64))
            basis = table[0]
            if sum(e.value_shape()) == 1:
                for i, value in enumerate(basis[0]):
                    assert np.isclose(value, reference[i](x))
            else:
                for i, ref in enumerate(reference):
                    assert np.allclose(basis[0][i::len(reference)], ref(x))
