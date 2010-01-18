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

from ufl import *
from ffc.finiteelement import FiniteElement as FFCElement
from ffc import jit

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

        P1 = FFCElement(FiniteElement("Lagrange", "triangle", 1))
        self.assertEqual(P1.space_dimension(), 3)

        P2 = FFCElement(FiniteElement("Lagrange", "triangle", 2))
        self.assertEqual(P2.space_dimension(), 6)

        P3 = FFCElement(FiniteElement("Lagrange", "triangle", 3))
        self.assertEqual(P3.space_dimension(), 10)

    def testDiscontinuousLagrange(self):
        "Test space dimensions of discontinuous Lagrange elements."

        P0 = FFCElement(FiniteElement("Discontinuous Lagrange", "triangle", 0))
        self.assertEqual(P0.space_dimension(), 1)

        P1 = FFCElement(FiniteElement("Discontinuous Lagrange", "triangle", 1))
        self.assertEqual(P1.space_dimension(), 3)

        P2 = FFCElement(FiniteElement("Discontinuous Lagrange", "triangle", 2))
        self.assertEqual(P2.space_dimension(), 6)

        P3 = FFCElement(FiniteElement("Discontinuous Lagrange", "triangle", 3))
        self.assertEqual(P3.space_dimension(), 10)

class FunctionValueTests(unittest.TestCase):

    def testContinuousLagrange1D(self):
        "Test values of continuous Lagrange functions in 1D."

        element = FFCElement(FiniteElement("Lagrange", "interval", 1))
        basis = element.basis()

        reference = [lambda x: 1 - x[0],
                     lambda x: x[0]]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(interval)
                self.assertAlmostEqual(basis[i](x), reference[i](x))

    def testContinuousLagrange2D(self):
        "Test values of continuous Lagrange functions in 2D."

        element = FFCElement(FiniteElement("Lagrange", "triangle", 1))
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

        element = FFCElement(FiniteElement("Lagrange", "tetrahedron", 1))
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

        element = FFCElement(FiniteElement("Discontinuous Lagrange", "interval", 1))
        basis = element.basis()

        reference = [lambda x: 1 - x[0],
                     lambda x: x[0]]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(interval)
                self.assertAlmostEqual(basis[i](x), reference[i](x))

    def testDiscontinuousLagrange2D(self):
        "Test values of discontinuous Lagrange functions in 2D."

        element = FFCElement(FiniteElement("Discontinuous Lagrange", "triangle", 1))
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

        element = FFCElement(FiniteElement("Discontinuous Lagrange", "tetrahedron", 1))
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

        element = FFCElement(FiniteElement("Brezzi-Douglas-Marini", "triangle", 1))
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

    def testRT1(self):
        "Test values of RT1."

        element = FFCElement(FiniteElement("Raviart-Thomas", "triangle", 1))
        basis = element.basis()

        reference = [lambda x: (x[0], x[1]),
                     lambda x: (1 - x[0], -x[1]),
                     lambda x: (x[0], x[1] - 1)]

        for i in range(len(basis)):
            for j in range(num_points):
                x = random_point(triangle)
                self.assertAlmostEqual(basis[i](x)[0], reference[i](x)[0])
                self.assertAlmostEqual(basis[i](x)[1], reference[i](x)[1])

    def testRT2(self):
        "Test values of RT2."

        element = FFCElement(FiniteElement("Raviart-Thomas", "triangle", 2))
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

        # FIXME: Use local cache: cache_dir argument to instant.build_module
        options = {"log_level": INFO + 5}

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
    os.system("python testcreateelement.py")
    os.system("python %s" % os.path.join("symbolics", "testsymbolics.py"))
    unittest.main()
