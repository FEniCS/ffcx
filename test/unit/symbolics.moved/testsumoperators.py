#!/usr/bin/env python

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-01-06"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-07

# Pyhton modules
import unittest
import time

# FFC modules
from ffc.quadrature.reduce_operations import operation_count, expand_operations, reduce_operations
from ffc.quadrature.symbolics import *
from ffc.quadrature.sumobj import _group_fractions
from ffc.ufcformat import Format
from ffc.constants import FFC_OPTIONS
from ffc.log import error, push_level, pop_level, CRITICAL

class TestSumOperators(unittest.TestCase):

    def testSumOperators(self):
        "Test binary operators"

        f2 = FloatValue(2.0)
        fm3 = FloatValue(-3.0)

        x = Symbol("x", GEO)
        y = Symbol("y", GEO)
        z = Symbol("z", GEO)

        p0 = Product([f2, x])
        p1 = Product([x, y])

        S0 = Sum([x, y])
        S1 = Sum([x, z])

        F0 = Fraction(p0, y)

        # Test Sum '+'
        self.assertEqual(str(S0 + f2), '(2 + x + y)')
        self.assertEqual(str(S0 + x), '(x + x + y)')
        self.assertEqual(str(S0 + p0), '(x + y + 2*x)')
        self.assertEqual(str(S0 + S0), '(x + x + y + y)')
        self.assertEqual(str(S0 + F0), '(x + y + 2*x/y)')

        # Test Sum '-'
        self.assertEqual(str(S0 - f2), '(x + y-2)')
        self.assertEqual(str(S0 - fm3), '(x + y + 3)')
        self.assertEqual(str(S0 - x), '(x + y - x)')
        self.assertEqual(str(S0 - p0), '(x + y-2*x)')
        self.assertEqual(str(S0 - Product([fm3, p0])), '(x + y + 6*x)')
        self.assertEqual(str(S0 - S0), '(x + y - (x + y))')
        self.assertEqual(str(S0 - F0), '(x + y - 2*x/y)')

        # Test Sum '*'
        self.assertEqual(str(S0 * f2), '(2*x + 2*y)')
        self.assertEqual(str(S0 * x), '(x*x + x*y)')
        self.assertEqual(str(S0 * p0), '(2*x*x + 2*x*y)')
        self.assertEqual(str(S0 * S0), '(2*x*y + x*x + y*y)')
        self.assertEqual(str(S0 * F0), '(2*x + 2*x*x/y)')

        # Test Sum '/'
        self.assertEqual(str(S0 / f2), '(0.5*x + 0.5*y)')
        self.assertEqual(str(S0 / x), '(1 + y/x)')
        self.assertEqual(str(S0 / p0), '(0.5 + 0.5*y/x)')
        self.assertEqual(str(S0 / p1), '(1/x + 1/y)')
        self.assertEqual(str(S0 / S0), '(x + y)/(x + y)')
        self.assertEqual(str(S0 / S1), '(x + y)/(x + z)')
        # Silence output
        push_level(CRITICAL)
        self.assertRaises(Exception, S0.__div__, FloatValue(0))
        self.assertRaises(Exception, S0.__div__, F0)
        pop_level()

if __name__ == "__main__":

    if format == None:
        set_format(Format(FFC_OPTIONS).format)

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestSumOperators('testSumOperators'))

