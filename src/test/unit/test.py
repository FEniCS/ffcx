"Unit tests for FFC"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-06 -- 2007-10-16"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

import unittest
import sys

sys.path.append("../..")
from ffc import *

class FiniteElementTests(unittest.TestCase):

    def testContinuousLagrange(self):
        "Test creation of continuous Lagrange elements"
        
        P1 = FiniteElement("Lagrange", "triangle", 1)
        self.assertEqual(P1.space_dimension(), 3)

        P2 = FiniteElement("Lagrange", "triangle", 2)
        self.assertEqual(P2.space_dimension(), 6)
        
        P3 = FiniteElement("Lagrange", "triangle", 3)
        self.assertEqual(P3.space_dimension(), 10)

    def testDiscontinuousLagrange(self):
        "Test creation of discontinuous Lagrange elements"

        P0 = FiniteElement("Discontinuous Lagrange", "triangle", 0)
        self.assertEqual(P0.space_dimension(), 1)
        
        P1 = FiniteElement("Discontinuous Lagrange", "triangle", 1)
        self.assertEqual(P1.space_dimension(), 3)

        P2 = FiniteElement("Discontinuous Lagrange", "triangle", 2)
        self.assertEqual(P2.space_dimension(), 6)
        
        P3 = FiniteElement("Discontinuous Lagrange", "triangle", 3)
        self.assertEqual(P3.space_dimension(), 10)

    def testFunctionValues2D(self):
        "Test values of simple Lagrange functions in 2D"

        P1 = FiniteElement("Lagrange", "triangle", 1)
        basis = P1.basis()

        x = [(-1, -1), (1, -1), (-1, 1)]
        for i in range(3):
            for j in range(3):
                value = basis[i](x[j])
                if i == j:
                    self.assertAlmostEqual(value, 1)
                else:
                    self.assertAlmostEqual(value, 0)
                    
    def testFunctionValues3D(self):
        "Test values of simple Lagrange functions in 3D"

        P1 = FiniteElement("Lagrange", "tetrahedron", 1)
        basis = P1.basis()

        x = [(-1, -1, -1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
        for i in range(4):
            for j in range(4):
                value = basis[i](x[j])
                if i == j:
                    self.assertAlmostEqual(value, 1)
                else:
                    self.assertAlmostEqual(value, 0)

if __name__ == "__main__":
    unittest.main()
