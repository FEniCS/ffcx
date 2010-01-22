#!/usr/bin/env python

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-01-06"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-06

# Pyhton modules
import unittest
import time

# FFC modules
from ffc.quadrature.symbolics import *
from ffc.cpp import format

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

    if format == None:
        set_format(format)

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestReduceGIP('testReduceGIP'))

