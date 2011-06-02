"Unit tests for FFC"

# Copyright (C) 2007-2009 Anders Logg
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
#
# First added:  2007-02-06
# Last changed: 2009-02-24

import unittest
import sys
import numpy
import math
import os
import instant
from time import time

sys.path.append(os.path.join(os.pardir, os.pardir))

from ufl import *
from ffc.fiatinterface import create_element as create
from ffc import jit

interval = [(0,), (1,)]
triangle = [(0, 0), (1, 0), (0, 1)]
tetrahedron = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
num_points = 5

# FIXME: A hack that will make the test run on Ubuntu 4.11, which ships two
# FIXME: versions of SWIG. UFC automatically uses the swig2.0 binary, which
# FIXME: FFC also need to use
use_swig2_binary = instant.check_and_set_swig_binary("swig2.0")

def random_point(shape):
    w = numpy.random.random(len(shape))
    return sum([numpy.array(shape[i])*w[i] for i in range(len(shape))]) / sum(w)

class SpaceDimensionTests(unittest.TestCase):

    def testContinuousLagrange(self):
        "Test space dimensions of continuous Lagrange elements."

        P1 = create(FiniteElement("Lagrange", "triangle", 1))
        self.assertEqual(P1.space_dimension(), 3)

        P2 = create(FiniteElement("Lagrange", "triangle", 2))
        self.assertEqual(P2.space_dimension(), 6)

        P3 = create(FiniteElement("Lagrange", "triangle", 3))
        self.assertEqual(P3.space_dimension(), 10)

    def testDiscontinuousLagrange(self):
        "Test space dimensions of discontinuous Lagrange elements."

        P0 = create(FiniteElement("DG", "triangle", 0))
        self.assertEqual(P0.space_dimension(), 1)

        P1 = create(FiniteElement("DG", "triangle", 1))
        self.assertEqual(P1.space_dimension(), 3)

        P2 = create(FiniteElement("DG", "triangle", 2))
        self.assertEqual(P2.space_dimension(), 6)

        P3 = create(FiniteElement("DG", "triangle", 3))
        self.assertEqual(P3.space_dimension(), 10)

class FunctionValueTests(unittest.TestCase):
    """
    These tests examine tabulate gives the correct answers for a the
    supported (non-mixed) elements of polynomial degree less than or
    equal to 3
    """

    # FIXME: Add tests for NED and BDM/RT in 3D.

    def _check_function_values(self, points, element, reference):
        for x in points:
            table = element.tabulate(0, (x,))
            basis = table[table.keys()[0]]
            for i in range(len(basis)):
                if element.value_shape() == ():
                    self.assertAlmostEqual(basis[i], reference[i](x))
                else:
                    for k in range(element.value_shape()[0]):
                        self.assertAlmostEqual(basis[i][k][0],
                                               reference[i](x)[k])

    def testContinuousLagrange1D(self):
        "Test values of continuous Lagrange functions in 1D."

        element = create(FiniteElement("Lagrange", "interval", 1))
        reference = [lambda x: 1 - x[0],
                     lambda x: x[0]]

        points = [random_point(interval) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testContinuousLagrange2D(self):
        "Test values of continuous Lagrange functions in 2D."

        element = create(FiniteElement("Lagrange", "triangle", 1))
        reference = [lambda x: 1 - x[0] - x[1],
                     lambda x: x[0],
                     lambda x: x[1]]

        points = [random_point(triangle) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testContinuousLagrange3D(self):
        "Test values of continuous Lagrange functions in 3D."

        element = create(FiniteElement("Lagrange", "tetrahedron", 1))
        reference = [lambda x: 1 - x[0] - x[1] - x[2],
                     lambda x: x[0],
                     lambda x: x[1],
                     lambda x: x[2]]

        points = [random_point(tetrahedron) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testDiscontinuousLagrange1D(self):
        "Test values of discontinuous Lagrange functions in 1D."

        element = create(FiniteElement("DG", "interval", 1))
        reference = [lambda x: 1 - x[0],
                     lambda x: x[0]]

        points = [random_point(interval) for i in range(num_points)]
        self._check_function_values(points, element, reference)


    def testDiscontinuousLagrange2D(self):
        "Test values of discontinuous Lagrange functions in 2D."

        element = create(FiniteElement("DG", "triangle", 1))
        reference = [lambda x: 1 - x[0] - x[1],
                     lambda x: x[0],
                     lambda x: x[1]]

        points = [random_point(triangle) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testDiscontinuousLagrange3D(self):
        "Test values of discontinuous Lagrange functions in 3D."

        element = create(FiniteElement("DG", "tetrahedron", 1))
        reference = [lambda x: 1 - x[0] - x[1] - x[2],
                     lambda x: x[0],
                     lambda x: x[1],
                     lambda x: x[2]]

        points = [random_point(tetrahedron) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testBDM1_2D(self):
        "Test values of BDM1."

        element = create(FiniteElement("Brezzi-Douglas-Marini", "triangle", 1))
        reference = [lambda x: (2*x[0], -x[1]),
                     lambda x: (-x[0], 2*x[1]),
                     lambda x: (2 - 2*x[0] - 3*x[1], x[1]),
                     lambda x: (- 1 + x[0] + 3*x[1], - 2*x[1]),
                     lambda x: (-x[0], -2 + 3*x[0] + 2*x[1]),
                     lambda x: (2*x[0], 1 - 3*x[0] - x[1])]

        points = [random_point(triangle) for i in range(num_points)]
        self._check_function_values(points, element, reference)


    def testRT1_2D(self):
        "Test values of RT1."

        element = create(FiniteElement("Raviart-Thomas", "triangle", 1))
        reference = [lambda x: (x[0], x[1]),
                     lambda x: (1 - x[0], -x[1]),
                     lambda x: (x[0], x[1] - 1)]

        points = [random_point(triangle) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testRT2_2D(self):
        "Test values of RT2."

        element = create(FiniteElement("Raviart-Thomas", "triangle", 2))

        reference = [ lambda x: (-x[0] + 3*x[0]**2,
                                 -x[1] + 3*x[0]*x[1]),
                      lambda x: (-x[0] + 3*x[0]*x[1],
                                 -x[1] + 3*x[1]**2),
                      lambda x: ( 2 - 5*x[0] - 3*x[1] + 3*x[0]*x[1] + 3*x[0]**2,
                                  -2*x[1] + 3*x[0]*x[1] + 3*x[1]**2),
                      lambda x: (-1.0 + x[0] + 3*x[1] - 3*x[0]*x[1],
                                 x[1] - 3*x[1]**2),
                      lambda x: (2*x[0] - 3*x[0]*x[1] - 3*x[0]**2,
                                 -2 + 3*x[0]+ 5*x[1] - 3*x[0]*x[1] - 3*x[1]**2),
                      lambda x: (- x[0] + 3*x[0]**2,
                                 + 1 - 3*x[0] - x[1] + 3*x[0]*x[1]),
                      lambda x: (6*x[0] - 3*x[0]*x[1] - 6*x[0]**2,
                                 3*x[1] - 6*x[0]*x[1] - 3*x[1]**2),
                      lambda x: (3*x[0] - 6*x[0]*x[1] - 3*x[0]**2,
                                 6*x[1]- 3*x[0]*x[1] - 6*x[1]**2),
                      ]


        points = [random_point(triangle) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testNED1_2D(self):
        "Test values of NED1."

        element = create(FiniteElement("N1curl", "triangle", 1))
        reference = [ lambda x: (-x[1], x[0]),
                      lambda x: ( x[1], 1 - x[0]),
                      lambda x: ( 1 - x[1], x[0]),
                      ]

        points = [random_point(triangle) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testRT1_3D(self):
        element = create(FiniteElement("RT", "tetrahedron", 1))
        reference = [lambda x: (-x[0], -x[1], -x[2]),
                     lambda x: (-1.0 + x[0], x[1], x[2]),
                     lambda x: (-x[0], 1.0 - x[1], -x[2]),
                     lambda x: ( x[0], x[1], -1.0 + x[2])
                     ]
        points = [random_point(tetrahedron) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testBDM1_3D(self):
        element = create(FiniteElement("BDM", "tetrahedron", 1))
        reference = [ lambda x: (-3*x[0], x[1], x[2]),
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
                      lambda x: (-x[0], 3*x[1], 1.0 - 4*x[1] - x[2])
                      ]
        points = [random_point(tetrahedron) for i in range(num_points)]
        self._check_function_values(points, element, reference)

    def testNED1_3D(self):
        element = create(FiniteElement("N1curl", "tetrahedron", 1))
        reference = [ lambda x: (0.0, -x[2], x[1]),
                      lambda x: (-x[2], 0.0,  x[0]),
                      lambda x: (-x[1],  x[0], 0.0),
                      lambda x: ( x[2], x[2], 1.0 - x[0] - x[1]),
                      lambda x: (x[1], 1.0 - x[0] - x[2], x[1]),
                      lambda x: (1.0 - x[1] - x[2], x[0], x[0])
                      ]
        points = [random_point(tetrahedron) for i in range(num_points)]
        self._check_function_values(points, element, reference)

class JITTests(unittest.TestCase):

    def testPoisson(self):
        "Test that JIT compiler is fast enough."

        # FIXME: Use local cache: cache_dir argument to instant.build_module
        #options = {"log_level": INFO + 5}
        #options = {"log_level": 5}
        options = {"log_level": WARNING}

        # FIXME: A hack to get the unit test pass on Ubuntu 11.04, see above.
        if use_swig2_binary:
            options["swig_binary"] = "swig2.0"

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

        print ""

        # Compile a1 (should be fairly fast, using disk cache)
        t = time()
        jit(a1, options)
        dt1 = time() - t

        # Good values
        dt0_good = 0.005
        dt1_good = 0.01

        print ""
        print "JIT in-memory cache:", dt0
        print "JIT disk cache:     ", dt1
        print "Reasonable values are %g and %g" % (dt0_good, dt1_good)

        # Check times
        self.assertTrue(dt0 < 10*dt0_good)
        self.assertTrue(dt1 < 10*dt1_good)

if __name__ == "__main__":
    unittest.main()
