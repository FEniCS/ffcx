#!/usr/bin/env python

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-01-06"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-28

# Pyhton modules
import unittest
import time

# FFC modules
from ffc.quadrature.reduce_operations import operation_count, expand_operations, reduce_operations
from ffc.quadrature.symbolics import *
from ffc.quadrature.sumobj import _group_fractions
from ffc.cpp import format, set_float_formatting
from ffc.constants import FFC_OPTIONS
set_float_formatting(FFC_OPTIONS)
from ffc.log import error, push_level, pop_level, CRITICAL

class TestSumOperators(unittest.TestCase):

    def testSumOperators(self):
        "Test binary operators"

        f_0_5 = format["float"](0.5)
        f_1 = format["float"](1)
        f_2 = format["float"](2)
        f_3 = format["float"](3)
        f_6 = format["float"](6)
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
        self.assertEqual(str(S0 + f2), '(%s + x + y)' % f_2)
        self.assertEqual(str(S0 + x), '(x + x + y)')
        self.assertEqual(str(S0 + p0), '(x + y + %s*x)' % f_2)
        self.assertEqual(str(S0 + S0), '(x + x + y + y)')
        self.assertEqual(str(S0 + F0), '(x + y + %s*x/y)' % f_2)

        # Test Sum '-'
        self.assertEqual(str(S0 - f2), '(x + y-%s)' % f_2)
        self.assertEqual(str(S0 - fm3), '(x + y + %s)' % f_3)
        self.assertEqual(str(S0 - x), '(x + y - x)')
        self.assertEqual(str(S0 - p0), '(x + y-%s*x)' % f_2)
        self.assertEqual(str(S0 - Product([fm3, p0])), '(x + y + %s*x)' % f_6)
        self.assertEqual(str(S0 - S0), '(x + y - (x + y))')
        self.assertEqual(str(S0 - F0), '(x + y - %s*x/y)' % f_2)

        # Test Sum '*'
        self.assertEqual(str(S0 * f2), '(%s*x + %s*y)' % (f_2, f_2))
        self.assertEqual(str(S0 * x), '(x*x + x*y)')
        self.assertEqual(str(S0 * p0), '(%s*x*x + %s*x*y)' % (f_2, f_2))
        self.assertEqual(str(S0 * S0), '(%s*x*y + x*x + y*y)' % f_2)
        self.assertEqual(str(S0 * F0), '(%s*x + %s*x*x/y)' % (f_2, f_2))

        # Test Sum '/'
        self.assertEqual(str(S0 / f2), '(%s*x + %s*y)' % (f_0_5, f_0_5))
        self.assertEqual(str(S0 / x), '(%s + y/x)' % f_1)
        self.assertEqual(str(S0 / p0), '(%s + %s*y/x)' % (f_0_5, f_0_5))
        self.assertEqual(str(S0 / p1), '(%s/x + %s/y)' % (f_1, f_1))
        self.assertEqual(str(S0 / S0), '(x + y)/(x + y)')
        self.assertEqual(str(S0 / S1), '(x + y)/(x + z)')
        # Silence output
        push_level(CRITICAL)
        self.assertRaises(Exception, S0.__div__, FloatValue(0))
        self.assertRaises(Exception, S0.__div__, F0)
        pop_level()

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestSumOperators('testSumOperators'))

