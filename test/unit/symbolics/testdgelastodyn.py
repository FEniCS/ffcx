#!/usr/bin/env python

# Copyright (C) 2010 Kristian B. Oelgaard
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2010-01-06
# Last changed: 2010-02-01

# Pyhton modules
import unittest
import time

# FFC modules
from ffc.quadrature.symbolics import *
from ffc.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])

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

