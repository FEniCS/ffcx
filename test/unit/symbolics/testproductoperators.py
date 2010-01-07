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

class TestProductOperators(unittest.TestCase):

    def testProductOperators(self):
        "Test binary operators"

        f0 = FloatValue(0.0)
        f1 = FloatValue(1.0)
        f2 = FloatValue(2.0)
        fm1 = FloatValue(-1.0)
        fm3 = FloatValue(-3.0)

        x = Symbol("x", GEO)
        y = Symbol("y", GEO)
        z = Symbol("z", GEO)

        p0 = Product([f2, x])
        p1 = Product([x, y])
        p2 = Product([f2, z])
        p3 = Product([x, y, z])

        S0 = Sum([x, y])
        S1 = Sum([x, z])

        F0 = Fraction(f2, x)
        F1 = Fraction(x, y)
        F2 = Fraction(x, S0)
        F3 = Fraction(x, y)
        F4 = Fraction(p0, y)

        # Test Product '+'
        self.assertEqual(str(p0 + f2), '(2 + 2*x)')
        self.assertEqual(str(p0 + x), '3*x')
        self.assertEqual(str(p0 + y), '(y + 2*x)')
        self.assertEqual(str(p0 + p0), '4*x')
        self.assertEqual(str(p0 + p1), '(2*x + x*y)')
        self.assertEqual(p0 + Product([fm1, x]), x)
        self.assertEqual(Product([fm1, x]) + x, f0)
        self.assertEqual(str(x + Product([fm1, x])), '0')
        self.assertEqual(str(p0 + S0), '(x + y + 2*x)')
        self.assertEqual(str(p0 + F3), '(2*x + x/y)')

        # Test Product '-'
        self.assertEqual(str(p0 - f2), '(2*x-2)')
        self.assertEqual(str(p0 - x), 'x')
        self.assertEqual(str(p0 - y), '(2*x - y)')
        self.assertEqual(str(p0 - p0), '0')
        self.assertEqual(str(p0 - p1), '(2*x - x*y)')
        self.assertEqual(str(p0 - S0), '(2*x - (x + y))')
        self.assertEqual(str(p0 - F3), '(2*x - x/y)')

        # Test Product '*', only need to test float, symbol and product.
        # Sum and fraction are handled by 'other'
        self.assertEqual(str(p0 * f0), '0')
        self.assertEqual(str(p0 * fm3), '-6*x')
        self.assertEqual(str(p0 * y), '2*x*y')
        self.assertEqual(str(p0 * p1), '2*x*x*y')

        # Test Product '/'
        self.assertEqual(str(Product([f0, x])/x), '0')
        self.assertEqual(str(p0/S0), '2*x/(x + y)')
        self.assertEqual(p1/y, x)
        self.assertEqual(p1/p2, Fraction(Product([p1, FloatValue(0.5)]), z))
        self.assertEqual(p1/z, Fraction(p1, z))
        self.assertEqual(p0/Product([f2, p1]), Fraction(f1, y))
        self.assertEqual(p1/p0,
        Product([FloatValue(0.5), y]))
        self.assertEqual(p1/p1, f1)
        self.assertEqual(p1/p3, Fraction(f1, z))
        self.assertEqual(str(p1/p3), '1/z')
        # Silence output
        push_level(CRITICAL)
        self.assertRaises(Exception, p0.__div__, f0)
        self.assertRaises(Exception, p0.__div__, F0)
        pop_level()

if __name__ == "__main__":

    if format == None:
        set_format(Format(FFC_OPTIONS).format)

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestProductOperators('testProductOperators'))

