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
from ffc.quadrature.symbolics import *
from ffc.cpp import format, set_float_formatting
from ffc.constants import FFC_OPTIONS
set_float_formatting(FFC_OPTIONS['precision'])

class TestDGElastoDyn(unittest.TestCase):

    def testDGElastoDyn(self):
        expr = Product([
                       Sum([
                            Symbol("F0", IP),
                            Symbol("F1", IP)
                          ]),
                       Fraction(
                                 Symbol("w4", GEO),
                                 Symbol("w3", GEO)
                                ),
                       Fraction(
                                 Product([
                                          Symbol("w2", GEO),
                                          Symbol("w5", GEO)
                                         ]),
                                 Symbol("w6", GEO)
                                )
                      ])

#        print "\nDGElastoDyn"
#        start = time.time()
        expr_exp = expr.expand()
#        print "DGElastoDyn: time, expand():     ", time.time() - start

#        start = time.time()
        expr_red = expr_exp.reduce_ops()
#        print "DGElastoDyn: time, reduce_ops(): ", time.time() - start

#        print "expr.ops():     ", expr.ops()
#        print "expr_exp.ops(): ", expr_exp.ops()
#        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        F0, F1, w2, w3, w4, w5, w6 = (3.12, -8.1, -45.3, 17.5, 2.2, 5.3, 9.145)
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertEqual(expr.ops(), 6)
        self.assertEqual(expr_exp.ops(), 11)
        self.assertEqual(expr_red.ops(), 6)

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestDGElastoDyn('testDGElastoDyn'))

