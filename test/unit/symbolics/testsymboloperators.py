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

class TestSymbolOperators(unittest.TestCase):

    def testSymbolOperators(self):
        "Test binary operators"

        f0 = FloatValue(0.0)
        f2 = FloatValue(2.0)
        fm1 = FloatValue(-1.0)
        fm3 = FloatValue(-3.0)

        x = Symbol("x", GEO)
        y = Symbol("y", GEO)
        z = Symbol("z", GEO)

        p0 = Product([f2, x])
        p1 = Product([x, y])
        p2 = Product([f2, z])
        p3 = Product([y, x, z])

        S0 = Sum([x, y])
        S1 = Sum([x, z])

        F0 = Fraction(f2, y)
        F1 = Fraction(x, y)
        F2 = Fraction(x, S0)
        F3 = Fraction(x, y)
        F4 = Fraction(p0, y)
        F5 = Fraction(fm3, y)

        # Test Symbol '+'
        self.assertEqual(str(x + f2), '(2 + x)')
        self.assertEqual(str(x + x), '2*x')
        self.assertEqual(str(x + y), '(x + y)')
        self.assertEqual(str(x + p0), '3*x')
        self.assertEqual(str(x + p1), '(x + x*y)')
        self.assertEqual(str(x + S0), '(x + x + y)')
        self.assertEqual(str(x + F0), '(x + 2/y)')

        # Test Symbol '-'
        self.assertEqual(str(x - f2), '(x-2)')
        self.assertEqual(str(x - x), '0')
        self.assertEqual(str(x - y), '(x - y)')
        self.assertEqual(str(x - p0), ' - x')
        self.assertEqual(str(x - p1), '(x - x*y)')
        self.assertEqual(str(x - S0), '(x - (x + y))')
        self.assertEqual(str(x - F5), '(x - -3/y)')

        # Test Symbol '*', only need to test float, symbol and product. Sum and
        # fraction are handled by 'other'
        self.assertEqual(str(x*f2), '2*x')
        self.assertEqual(str(x*y), 'x*y')
        self.assertEqual(str(x*p1), 'x*x*y')

        # Test Symbol '/'
        self.assertEqual(str(x/f2), '0.5*x')
        self.assertEqual(str(x/x), '1')
        self.assertEqual(str(x/y), 'x/y')
        self.assertEqual(str(x/S0), 'x/(x + y)')
        self.assertEqual(str(x/p0), '0.5')
        self.assertEqual(str(y/p1), '1/x')
        self.assertEqual(str(z/p0), '0.5*z/x')
        self.assertEqual(str(z/p1), 'z/(x*y)')
        # Silence output
        push_level(CRITICAL)
        self.assertRaises(Exception, x.__div__, F0)
        self.assertRaises(Exception, y.__div__, FloatValue(0))
        pop_level()

if __name__ == "__main__":

    if format == None:
        set_format(Format(FFC_OPTIONS).format)

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestSymbolOperators('testSymbolOperators'))

