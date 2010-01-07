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

class TestFractionOperators(unittest.TestCase):

    def testFractionOperators(self):
        "Test binary operators"

        f2 = FloatValue(2.0)
        fm3 = FloatValue(-3.0)

        x = Symbol("x", GEO)
        y = Symbol("y", GEO)

        p0 = Product([f2, x])
        p1 = Product([x, y])

        S0 = Sum([x, y])

        F0 = Fraction(f2, y)
        F1 = Fraction(x, y)
        F2 = Fraction(x, S0)
        F3 = Fraction(x, y)
        F4 = Fraction(p0, y)
        F5 = Fraction(Product([fm3, x]), y)

        # Test Fraction '+'
        self.assertEqual(str(F0 + f2), '(2 + 2/y)')
        self.assertEqual(str(F1 + x), '(x + x/y)')
        self.assertEqual(str(F1 + p0), '(2*x + x/y)')
        self.assertEqual(str(F1 + S0), '(x + y + x/y)')
        self.assertEqual(str(F1 + F3), '2*x/y')
        self.assertEqual(str(F0 + F1), '(2 + x)/y')
        self.assertEqual(str(F2 + F4), '(2*x/y + x/(x + y))')

        # Test Fraction '-'
        self.assertEqual(str(F0 - f2), '(2/y-2)')
        self.assertEqual(str(F1 - x), '(x/y - x)')
        self.assertEqual(str(F1 - p0), '(x/y-2*x)')
        self.assertEqual(str(F1 - S0), '(x/y - (x + y))')
        self.assertEqual(str(F1 - F3), '0')
        self.assertEqual(str(F4 - F1), 'x/y')
        self.assertEqual(str(F4 - F5), '5*x/y')
        self.assertEqual(str(F0 - F1), '(2 - x)/y')
        self.assertEqual(str(F2 - F4), '(x/(x + y) - 2*x/y)')

        # Test Fraction '*'
        self.assertEqual(str(F1 * f2), '2*x/y')
        self.assertEqual(str(F1 * x), 'x*x/y')
        self.assertEqual(str(F1 * p1), 'x*x')
        self.assertEqual(str(F1 * S0), '(x + x*x/y)')
        self.assertEqual(repr(F1 * S0), repr(Sum([x, Fraction( Product([x, x]), y)]) ))
        self.assertEqual(str(F1 * F0), '2*x/(y*y)')

        # Test Fraction '/'
        self.assertEqual(str(F0 / f2), '1/y')
        self.assertEqual(str(F1 / x), '1/y')
        self.assertEqual(str(F4 / p1), '2/(y*y)')
        self.assertEqual(str(F4 / x), '2/y')
        self.assertEqual(str(F2 / y), 'x/(x*y + y*y)')
        self.assertEqual(str(F0 / S0), '2/(x*y + y*y)')
        # Silence output
        push_level(CRITICAL)
        self.assertRaises(Exception, F0.__div__, F0)
        pop_level()

if __name__ == "__main__":

    if format == None:
        set_format(Format(FFC_OPTIONS).format)

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestFractionOperators('testFractionOperators'))

