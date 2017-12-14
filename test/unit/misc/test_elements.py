# -*- coding: utf-8 -*-
"Unit tests for FFC"

# Copyright (C) 2007-2017 Anders Logg and Garth N. Wells
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Marie E. Rognes, 2010
# Modified by Lizao Li, 2016


import pytest
import os
import sys
import numpy
import math
from time import time

from ufl import *
from ffc.fiatinterface import create_element
from ffc import jit


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
    w = numpy.random.random(len(shape))
    return sum([numpy.array(shape[i])*w[i] for i in range(len(shape))])/sum(w)


@pytest.mark.parametrize("degree, expected_dim", [(1, 3), (2, 6), (3, 10)])
def test_continuous_lagrange(degree, expected_dim):
    "Test space dimensions of continuous Lagrange elements."
    P = create_element(FiniteElement("Lagrange", "triangle", degree))
    assert P.space_dimension() == expected_dim

@pytest.mark.parametrize("degree, expected_dim", [(1, 4), (2, 9), (3, 16)])
def test_continuous_lagrange_quadrilateral(degree, expected_dim):
    "Test space dimensions of continuous TensorProduct elements (quadrilateral)."
    P = create_element(FiniteElement("Lagrange", "quadrilateral", degree))
    assert P.space_dimension() == expected_dim

@pytest.mark.parametrize("degree, expected_dim", [(0, 1), (1, 3), (2, 6), (3, 10)])
def test_discontinuous_lagrange(degree, expected_dim):
    "Test space dimensions of discontinuous Lagrange elements."
    P = create_element(FiniteElement("DG", "triangle", degree))
    assert P.space_dimension() == expected_dim

@pytest.mark.parametrize("degree, expected_dim",
                         [(0, 3), (1, 9), (2, 18), (3, 30)])
def test_regge(degree, expected_dim):
    "Test space dimensions of generalized Regge element."
    P = create_element(FiniteElement("Regge", "triangle", degree))
    assert P.space_dimension() == expected_dim

@pytest.mark.parametrize("degree, expected_dim",
                         [(0, 3), (1, 9), (2, 18), (3, 30)])
def test_hhj(degree, expected_dim):
    "Test space dimensions of Hellan-Herrmann-Johnson element."
    P = create_element(FiniteElement("HHJ", "triangle", degree))
    assert P.space_dimension() == expected_dim


class TestFunctionValues():
    """These tests examine tabulate gives the correct answers for a the
supported (non-mixed) for low degrees"""

    # FIXME: Add tests for NED and BDM/RT in 3D.

    # Shape (basis) functions on reference element
    reference_interval_1 = [lambda x: 1 - x[0], lambda x: x[0]]
    reference_triangle_1 = [lambda x: 1 - x[0] - x[1], lambda x: x[0], lambda x: x[1]]
    reference_tetrahedron_1 = [lambda x: 1 - x[0] - x[1] - x[2], lambda x: x[0],
                               lambda x: x[1], lambda x: x[2]]
    reference_triangle_bdm1 = [lambda x: (2*x[0], -x[1]),
                               lambda x: (-x[0], 2*x[1]),
                               lambda x: (2 - 2*x[0] - 3*x[1], x[1]),
                               lambda x: (- 1 + x[0] + 3*x[1], - 2*x[1]),
                               lambda x: (-x[0], -2 + 3*x[0] + 2*x[1]),
                               lambda x: (2*x[0], 1 - 3*x[0] - x[1])]
    reference_triangle_rt1 = [lambda x: (x[0], x[1]), lambda x: (1 - x[0], -x[1]),
                              lambda x: (x[0], x[1] - 1)]
    reference_triangle_rt2 = [lambda x: (-x[0] + 3*x[0]**2, -x[1] + 3*x[0]*x[1]),
                              lambda x: (-x[0] + 3*x[0]*x[1], -x[1] + 3*x[1]**2),
                              lambda x: ( 2 - 5*x[0] - 3*x[1] + 3*x[0]*x[1] + 3*x[0]**2,
                                          -2*x[1] + 3*x[0]*x[1] + 3*x[1]**2),
                              lambda x: (-1.0 + x[0] + 3*x[1] - 3*x[0]*x[1], x[1] - 3*x[1]**2),
                              lambda x: (2*x[0] - 3*x[0]*x[1] - 3*x[0]**2,
                                         -2 + 3*x[0]+ 5*x[1] - 3*x[0]*x[1] - 3*x[1]**2),
                              lambda x: (- x[0] + 3*x[0]**2,
                                       + 1 - 3*x[0] - x[1] + 3*x[0]*x[1]),
                              lambda x: (6*x[0] - 3*x[0]*x[1] - 6*x[0]**2,
                                         3*x[1] - 6*x[0]*x[1] - 3*x[1]**2),
                              lambda x: (3*x[0] - 6*x[0]*x[1] - 3*x[0]**2,
                                         6*x[1]- 3*x[0]*x[1] - 6*x[1]**2)]
    reference_triangle_ned1 = [lambda x: (-x[1], x[0]), lambda x: ( x[1], 1 - x[0]),
                               lambda x: ( 1 - x[1], x[0])]
    reference_tetrahedron_rt1 = [lambda x: (-x[0], -x[1], -x[2]),
                                 lambda x: (-1.0 + x[0], x[1], x[2]),
                                 lambda x: (-x[0], 1.0 - x[1], -x[2]),
                                 lambda x: ( x[0], x[1], -1.0 + x[2])]
    reference_tetrahedron_bdm1 = [lambda x: (-3*x[0], x[1], x[2]),
                                  lambda x: (x[0], -3*x[1], x[2]),
                                  lambda x: (x[0], x[1], -3*x[2]),
                                  lambda x: (-3.0 + 3*x[0] + 4*x[1] + 4*x[2], -x[1], -x[2]),
                                  lambda x: (1.0 - x[0] - 4*x[1], 3*x[1], -x[2]),
                                  lambda x: (1.0 - x[0] - 4*x[2], -x[1], 3*x[2]),
                                  lambda x: (x[0], 3.0 - 4*x[0] - 3*x[1] - 4*x[2], x[2]),
                                  lambda x: (-3*x[0], -1.0 + 4*x[0] + x[1], x[2]),
                                  lambda x: (x[0], -1.0 + x[1] + 4*x[2], -3*x[2]),
                                  lambda x: (-x[0], -x[1], -3.0 + 4*x[0] + 4*x[1] + 3*x[2]),
                                  lambda x: (3*x[0], -x[1], 1.0 - 4*x[0] - x[2]),
                                  lambda x: (-x[0], 3*x[1], 1.0 - 4*x[1] - x[2])]
    reference_tetrahedron_ned1 = [lambda x: (0.0, -x[2], x[1]),
                                  lambda x: (-x[2], 0.0,  x[0]),
                                  lambda x: (-x[1],  x[0], 0.0),
                                  lambda x: ( x[2], x[2], 1.0 - x[0] - x[1]),
                                  lambda x: (x[1], 1.0 - x[0] - x[2], x[1]),
                                  lambda x: (1.0 - x[1] - x[2], x[0], x[0])]
    reference_quadrilateral_1 = [lambda x: (1-x[0])*(1-x[1]),
                                 lambda x: (1-x[0])*x[1],
                                 lambda x: x[0]*(1-x[1]),
                                 lambda x: x[0]*x[1]]
    reference_hexahedron_1 = [lambda x: (1-x[0])*(1-x[1])*(1-x[2]),
                              lambda x: (1-x[0])*(1-x[1])*x[2],
                              lambda x: (1-x[0])*x[1]*(1-x[2]),
                              lambda x: (1-x[0])*x[1]*x[2],
                              lambda x: x[0]*(1-x[1])*(1-x[2]),
                              lambda x: x[0]*(1-x[1])*x[2],
                              lambda x: x[0]*x[1]*(1-x[2]),
                              lambda x: x[0]*x[1]*x[2]]

    # Tests to perform
    tests = [("Lagrange", "interval", 1, reference_interval_1),
             ("Lagrange", "triangle", 1, reference_triangle_1),
             ("Lagrange", "tetrahedron", 1, reference_tetrahedron_1),
             ("Lagrange", "quadrilateral", 1, reference_quadrilateral_1),
             ("Lagrange", "hexahedron", 1, reference_hexahedron_1),
             ("Discontinuous Lagrange", "interval", 1, reference_interval_1),
             ("Discontinuous Lagrange", "triangle", 1, reference_triangle_1),
             ("Discontinuous Lagrange", "tetrahedron", 1, reference_tetrahedron_1),
             ("Brezzi-Douglas-Marini", "triangle", 1, reference_triangle_bdm1),
             ("Raviart-Thomas", "triangle", 1, reference_triangle_rt1),
             ("Raviart-Thomas", "triangle", 2, reference_triangle_rt2),
             ("Discontinuous Raviart-Thomas", "triangle", 1, reference_triangle_rt1),
             ("Discontinuous Raviart-Thomas", "triangle", 2, reference_triangle_rt2),
             ("N1curl", "triangle", 1, reference_triangle_ned1),
             ("Raviart-Thomas", "tetrahedron", 1, reference_tetrahedron_rt1),
             ("Discontinuous Raviart-Thomas", "tetrahedron", 1, reference_tetrahedron_rt1),
             ("Brezzi-Douglas-Marini", "tetrahedron", 1, reference_tetrahedron_bdm1),
             ("N1curl", "tetrahedron", 1, reference_tetrahedron_ned1),
        ]


    @pytest.mark.parametrize("family, cell, degree, reference", tests)
    def test_values(self, family, cell, degree, reference):
        # Create element
        element = create_element(FiniteElement(family, cell, degree))

        # Get some points and check basis function values at points
        points = [random_point(element_coords(cell)) for i in range(5)]
        for x in points:
            table = element.tabulate(0, (x,))
            basis = table[list(table.keys())[0]]
            for i in range(len(basis)):
                if not element.value_shape():
                    assert round(float(basis[i]) - reference[i](x), 10) == 0.0
                else:
                    for k in range(element.value_shape()[0]):
                        assert round(basis[i][k][0] - reference[i](x)[k] , 10) == 0.0


class Test_JIT():

    def test_poisson(self):
        "Test that JIT compiler is fast enough."

        # FIXME: Use local cache: cache_dir argument to instant.build_module
        #options = {"log_level": INFO + 5}
        #options = {"log_level": 5}
        options = {"log_level": WARNING}

        # Define two forms with the same signatures
        element = FiniteElement("Lagrange", "triangle", 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        f = Coefficient(element)
        g = Coefficient(element)
        a0 = f*dot(grad(v), grad(u))*dx
        a1 = g*dot(grad(v), grad(u))*dx

        # Strange this needs to be done twice

        # Compile a0 so it will be in the cache (both in-memory and disk)
        jit(a0, options)
        jit(a0, options)

        # Compile a0 again (should be really fast, using in-memory cache)
        t = time()
        jit(a0, options)
        dt0 = time() - t

        # Compile a1 (should be fairly fast, using disk cache)
        t = time()
        jit(a1, options)
        dt1 = time() - t

        # Good values
        dt0_good = 0.005
        dt1_good = 0.01

        print("\nJIT in-memory cache:", dt0)
        print("JIT disk cache:     ", dt1)
        print("Reasonable values are %g and %g" % (dt0_good, dt1_good))

        # Check times
        assert dt0 < 10*dt0_good
        assert dt1 < 10*dt1_good
