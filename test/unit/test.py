"Unit tests for FFC"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-06 -- 2009-02-24"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

import unittest
import sys
import numpy
import math
import os
from time import time

sys.path.append(os.path.join(os.pardir, os.pardir))
from ffc import *

interval = [(0,), (1,)]
triangle = [(0, 0), (1, 0), (0, 1)]
tetrahedron = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
num_points = 5

def random_point(shape):
    w = numpy.random.random(len(shape))
    return sum([numpy.array(shape[i])*w[i] for i in range(len(shape))]) / sum(w)

class SpaceDimensionTests(unittest.TestCase):

    def testContinuousLagrange(self):
        "Test space dimensions of continuous Lagrange elements."
        
        P1 = FiniteElement("Lagrange", "triangle", 1)
        self.assertEqual(P1.space_dimension(), 3)

        P2 = FiniteElement("Lagrange", "triangle", 2)
        self.assertEqual(P2.space_dimension(), 6)
        
        P3 = FiniteElement("Lagrange", "triangle", 3)
        self.assertEqual(P3.space_dimension(), 10)

    def testDiscontinuousLagrange(self):
        "Test space dimensions of discontinuous Lagrange elements."

        P0 = FiniteElement("Discontinuous Lagrange", "triangle", 0)
        self.assertEqual(P0.space_dimension(), 1)
        
        P1 = FiniteElement("Discontinuous Lagrange", "triangle", 1)
        self.assertEqual(P1.space_dimension(), 3)

        P2 = FiniteElement("Discontinuous Lagrange", "triangle", 2)
        self.assertEqual(P2.space_dimension(), 6)
        
        P3 = FiniteElement("Discontinuous Lagrange", "triangle", 3)
        self.assertEqual(P3.space_dimension(), 10)

class FunctionValueTests(unittest.TestCase):

    def testContinuousLagrange1D(self):
        "Test values of continuous Lagrange functions in 1D."

        element = FiniteElement("Lagrange", "interval", 1)
        basis = element.basis()
        
        reference = [lambda x: 1 - x[0],
                     lambda x: x[0]]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(interval)
                self.assertAlmostEqual(basis[i](x), reference[i](x))

    def testContinuousLagrange2D(self):
        "Test values of continuous Lagrange functions in 2D."

        element = FiniteElement("Lagrange", "triangle", 1)
        basis = element.basis()
        
        reference = [lambda x: 1 - x[0] - x[1],
                     lambda x: x[0],
                     lambda x: x[1]]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(triangle)
                self.assertAlmostEqual(basis[i](x), reference[i](x))

    def testContinuousLagrange3D(self):
        "Test values of continuous Lagrange functions in 3D."

        element = FiniteElement("Lagrange", "tetrahedron", 1)
        basis = element.basis()
        
        reference = [lambda x: 1 - x[0] - x[1] - x[2],
                     lambda x: x[0],
                     lambda x: x[1],
                     lambda x: x[2]]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(tetrahedron)
                self.assertAlmostEqual(basis[i](x), reference[i](x))

    def testDiscontinuousLagrange1D(self):
        "Test values of discontinuous Lagrange functions in 1D."

        element = FiniteElement("Discontinuous Lagrange", "interval", 1)
        basis = element.basis()
        
        reference = [lambda x: 1 - x[0],
                     lambda x: x[0]]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(interval)
                self.assertAlmostEqual(basis[i](x), reference[i](x))

    def testDiscontinuousLagrange2D(self):
        "Test values of discontinuous Lagrange functions in 2D."

        element = FiniteElement("Discontinuous Lagrange", "triangle", 1)
        basis = element.basis()
        
        reference = [lambda x: 1 - x[0] - x[1],
                     lambda x: x[0],
                     lambda x: x[1]]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(triangle)
                self.assertAlmostEqual(basis[i](x), reference[i](x))

    def testDiscontinuousLagrange3D(self):
        "Test values of discontinuous Lagrange functions in 3D."

        element = FiniteElement("Discontinuous Lagrange", "tetrahedron", 1)
        basis = element.basis()
        
        reference = [lambda x: 1 - x[0] - x[1] - x[2],
                     lambda x: x[0],
                     lambda x: x[1],
                     lambda x: x[2]]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(tetrahedron)
                self.assertAlmostEqual(basis[i](x), reference[i](x))

    def testBDM1(self):
        "Test values of BDM1."

        element = FiniteElement("Brezzi-Douglas-Marini", "triangle", 1)
        basis = element.basis()
        
        reference = [lambda x: (2*x[0], -x[1]),
                     lambda x: (-x[0], 2*x[1]),
                     lambda x: (2 - 2*x[0] - 3*x[1], x[1]),
                     lambda x: (- 1 + x[0] + 3*x[1], - 2*x[1]),
                     lambda x: (-x[0], -2 + 3*x[0] + 2*x[1]),
                     lambda x: (2*x[0], 1 - 3*x[0] - x[1])]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(triangle)
                self.assertAlmostEqual(basis[i](x)[0], reference[i](x)[0])
                self.assertAlmostEqual(basis[i](x)[1], reference[i](x)[1])

    def testRT0(self):
        "Test values of RT0."

        element = FiniteElement("Raviart-Thomas", "triangle", 0)
        basis = element.basis()
        
        reference = [lambda x: (x[0], x[1]),
                     lambda x: (1 - x[0], -x[1]),
                     lambda x: (x[0], x[1] - 1)]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(triangle)
                self.assertAlmostEqual(basis[i](x)[0], reference[i](x)[0])
                self.assertAlmostEqual(basis[i](x)[1], reference[i](x)[1])

    def testRT1(self):
        "Test values of RT1."

        element = FiniteElement("Raviart-Thomas", "triangle", 1)
        basis = element.basis()
        
        reference = [lambda x: (-2*x[0] + 4*x[0]**2, -x[1] + 4*x[0]*x[1]),
                     lambda x: (-x[0] + 4*x[0]*x[1], -2*x[1] + 4*x[1]**2),
                     lambda x: (2 - 6*x[0] - 3*x[1] + 4*x[0]*x[1] + 4*x[0]**2,
                                -3*x[1] + 4*x[0]*x[1] + 4*x[1]**2),
                     lambda x: (-1 + x[0] +3*x[1] -4*x[0]*x[1],
                                2*x[1] -4*x[1]**2),
                     lambda x: (3*x[0] -4*x[0]*x[1] - 4*x[0]**2,
                                -2 + 3*x[0] +6*x[1] - 4*x[0]*x[1] -4*x[1]**2),
                     lambda x: (-2*x[0] +4* x[0]**2,
                                1 -3*x[0] -x[1] + 4*x[0]*x[1]),
                     lambda x: (16/math.sqrt(2)*x[0] - 8/math.sqrt(2)*x[0]*x[1]
                                - 16/math.sqrt(2)*x[0]**2,
                                8/math.sqrt(2)*x[1] - 16/math.sqrt(2)*x[0]*x[1]
                                - 8/math.sqrt(2)*x[1]**2),
                     lambda x: (8/math.sqrt(2)*x[0] -16/math.sqrt(2)*x[0]*x[1]
                                - 8/math.sqrt(2)*x[0]**2,
                                16/math.sqrt(2)*x[1] - 8/math.sqrt(2)*x[0]*x[1]
                                - 16/math.sqrt(2)*x[1]**2)
                     ]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(triangle)
                self.assertAlmostEqual(basis[i](x)[0], reference[i](x)[0])
                self.assertAlmostEqual(basis[i](x)[1], reference[i](x)[1])

class JITTests(unittest.TestCase):

    def testPoisson(self):
        "Test that JIT compiler is fast enough."

        # Define two forms with the same signatures
        element = FiniteElement("Lagrange", "triangle", 1)
        v = TestFunction(element)
        u = TrialFunction(element)
        a0 = dot(grad(v), grad(u))*dx
        a1 = dot(grad(v), grad(u))*dx

        # Compile a0 once so it will be in the cache (both in-memory and disk)
        jit(a0)

        # Compile a0 again (should be really fast, using in-memory cache)
        t = time()
        jit(a0)
        dt0 = time() - t

        # Compile a1 (should be fairly, using disk cache)
        t = time()
        jit(a1)
        dt1 = time() - t

        print ""
        print "JIT in-memory cache:", dt0
        print "JIT disk cache:     ", dt1
        print "Reasonable values are 0.002 and 0.03"
        
        # Check times
        self.assertTrue(dt0 < 0.001)
        self.assertTrue(dt1 < 0.1)
        
if __name__ == "__main__":
    unittest.main()
