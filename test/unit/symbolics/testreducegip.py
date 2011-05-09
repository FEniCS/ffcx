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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC.  If not, see <http://www.gnu.org/licenses/>.
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

class TestReduceGIP(unittest.TestCase):

    def testReduceGIP(self):

        expr = Sum([
                    Product([
                              Symbol("F17", IP), Symbol("F17", IP), Symbol("F18", IP), Symbol("F18", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G0", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F11", IP), Symbol("F13", IP), Symbol("F13", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G1", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F11", IP), Symbol("F12", IP), Symbol("F13", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F11", IP), Symbol("F12", IP), Symbol("F20", IP),
                              Symbol("F3", IP), Symbol("F8", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F18", IP), Symbol("F20", IP), Symbol("F3", IP),
                              Symbol("F8", IP), Symbol("F9", IP), Symbol("W9", IP), Symbol("G3", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F17", IP), Symbol("F18", IP), Symbol("F20", IP),
                              Symbol("F8", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F20", IP), Symbol("F8", IP), Symbol("F8", IP),
                              Symbol("F9", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F13", IP), Symbol("F17", IP), Symbol("F18", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ]),
                    Product([
                              Symbol("F20", IP), Symbol("F8", IP), Symbol("F8", IP), Symbol("F9", IP),
                              Symbol("F9", IP), Symbol("W9", IP), Symbol("G5", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F12", IP), Symbol("F17", IP), Symbol("F18", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G6", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F10", IP), Symbol("F20", IP), Symbol("F3", IP),
                              Symbol("F8", IP), Symbol("F8", IP), Symbol("W9", IP), Symbol("G1", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F10", IP), Symbol("F20", IP), Symbol("F8", IP),
                              Symbol("F8", IP), Symbol("W9", IP), Symbol("G7", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F17", IP), Symbol("F18", IP), Symbol("F19", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F12", IP), Symbol("F17", IP), Symbol("F19", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F11", IP), Symbol("F13", IP), Symbol("F13", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G7", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F17", IP), Symbol("F19", IP), Symbol("F20", IP),
                              Symbol("F3", IP), Symbol("F8", IP), Symbol("W9", IP), Symbol("G8", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F17", IP), Symbol("F18", IP), Symbol("F18", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G5", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F17", IP), Symbol("F19", IP), Symbol("F20", IP),
                              Symbol("F8", IP), Symbol("W9", IP), Symbol("G9", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F17", IP), Symbol("F18", IP), Symbol("F20", IP),
                              Symbol("F3", IP), Symbol("F8", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F18", IP), Symbol("F20", IP), Symbol("F8", IP),
                              Symbol("F9", IP), Symbol("W9", IP), Symbol("G6", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F11", IP), Symbol("F13", IP), Symbol("F20", IP),
                              Symbol("F3", IP), Symbol("F8", IP), Symbol("W9", IP), Symbol("G8", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F11", IP), Symbol("F12", IP), Symbol("F13", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F11", IP), Symbol("F13", IP), Symbol("F20", IP),
                              Symbol("F8", IP), Symbol("W9", IP), Symbol("G9", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F13", IP), Symbol("F20", IP), Symbol("F3", IP),
                              Symbol("F8", IP), Symbol("F9", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F12", IP), Symbol("F17", IP), Symbol("F18", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G3", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F20", IP), Symbol("F3", IP), Symbol("F8", IP),
                              Symbol("F8", IP), Symbol("F9", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F19", IP), Symbol("F20", IP), Symbol("F8", IP),
                              Symbol("F9", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F13", IP), Symbol("F17", IP), Symbol("F19", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G9", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F13", IP), Symbol("F17", IP), Symbol("F18", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F11", IP), Symbol("F12", IP), Symbol("F12", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G0", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F17", IP), Symbol("F19", IP), Symbol("F19", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G7", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F17", IP), Symbol("F18", IP), Symbol("F19", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("F8", IP), Symbol("F8", IP),
                              Symbol("F9", IP), Symbol("F9", IP), Symbol("W9", IP), Symbol("G0", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F12", IP), Symbol("F20", IP), Symbol("F8", IP),
                              Symbol("F9", IP), Symbol("W9", IP), Symbol("G6", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F13", IP), Symbol("F17", IP), Symbol("F19", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G8", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F19", IP), Symbol("F20", IP), Symbol("F3", IP),
                              Symbol("F8", IP), Symbol("F9", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ]),
                    Product([
                              Symbol("F10", IP), Symbol("F11", IP), Symbol("F12", IP), Symbol("F20", IP),
                              Symbol("F8", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F13", IP), Symbol("F20", IP), Symbol("F8", IP),
                              Symbol("F9", IP), Symbol("W9", IP), Symbol("G4", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F11", IP), Symbol("F12", IP), Symbol("F12", IP),
                              Symbol("F20", IP), Symbol("W9", IP), Symbol("G5", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F12", IP), Symbol("F20", IP), Symbol("F3", IP),
                              Symbol("F8", IP), Symbol("F9", IP), Symbol("W9", IP), Symbol("G3", GEO)
                            ]),
                    Product([
                              Symbol("F17", IP), Symbol("F17", IP), Symbol("F19", IP), Symbol("F19", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G1", GEO)
                            ]),
                    Product([
                              Symbol("F11", IP), Symbol("F12", IP), Symbol("F17", IP), Symbol("F19", IP),
                              Symbol("F20", IP), Symbol("F3", IP), Symbol("W9", IP), Symbol("G2", GEO)
                            ])
                   ])

#        print "\nReduceGIP"
#        start = time.time()
        expr_exp = expr.expand()
#        print "ReduceGIP: time, expand()      ", time.time() - start

#        start = time.time()
        expr_red = expr_exp.reduce_ops()
#        print "ReduceGIP: time, reduce_ops(): ", time.time() - start

#        print "expr.ops():     ", expr.ops()
#        print "expr_exp.ops(): ", expr_exp.ops()
#        print "expr_red.ops(): ", expr_red.ops()

#        print "expr: ", expr
#        print "exp:  ", expr_exp
#        print "red:  ", expr_red

        W9 = 9
        F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20 = [0.123 * i for i in range(1,21)]
        G0, G1, G2, G3, G4, G5, G6, G7, G8, G9 = [2.64 + 1.0/i for i in range(20, 30)]

        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertEqual(expr.ops(), 314)
        self.assertEqual(expr_exp.ops(), 314)
        self.assertEqual(expr_red.ops(), 114)

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestReduceGIP('testReduceGIP'))

