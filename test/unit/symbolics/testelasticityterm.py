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

class TestElasticityTerm(unittest.TestCase):

    def testElasticityTerm(self):
        # expr:  0.25*W1*det*(FE0_C2_D001[0][j]*FE0_C2_D001[0][k]*Jinv_00*Jinv_21 + FE0_C2_D001[0][j]*FE0_C2_D001[0][k]*Jinv_00*Jinv_21)
        expr = Product([
                         FloatValue(0.25), Symbol('W1', GEO), Symbol('det', GEO),
                         Sum([Product([Symbol('FE0_C2_D001_0_j', BASIS), Symbol('FE0_C2_D001_0_k', BASIS),
                                       Symbol('Jinv_00', GEO), Symbol('Jinv_21', GEO)]),
                              Product([Symbol('FE0_C2_D001_0_j', BASIS), Symbol('FE0_C2_D001_0_k', BASIS),
                                  Symbol('Jinv_00', GEO), Symbol('Jinv_21', GEO)])
                             ])
                      ])

#        print "\nElasticityTerm"
#        start = time.time()
        expr_exp = expr.expand()
#        print "ElasticityTerm: time, expand():     ", time.time() - start

#        start = time.time()
        expr_red = expr_exp.reduce_ops()
#        print "ElasticityTerm: time, reduce_ops(): ", time.time() - start

#        print "expr.ops():     ", expr.ops()
#        print "expr_exp.ops(): ", expr_exp.ops()
#        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        det, W1, Jinv_00, Jinv_21, FE0_C2_D001_0_j, FE0_C2_D001_0_k = [0.123 + i for i in range(6)]
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertEqual(expr.ops(), 10)
        self.assertEqual(expr_exp.ops(), 6)
        self.assertEqual(expr_red.ops(), 6)

        # Generate code
        ip_consts = {}
        geo_consts = {}
        trans_set = set()

        start = time.time()
        opt_code = optimise_code(expr, ip_consts, geo_consts, trans_set)
#        print "ElasticityTerm, optimise_code(): ", time.time() - start

        G0 = eval(str(geo_consts.items()[0][0]))
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestElasticityTerm('testElasticityTerm'))

