"Unit tests for FFC"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-06 -- 2007-02-06"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

import unittest
import sys

sys.path.append("..")
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

if __name__ == "__main__":
    unittest.main()
