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
from ffc.cpp import format, set_float_formatting
from ffc.constants import FFC_OPTIONS
set_float_formatting(FFC_OPTIONS['precision'])

class TestPoisson(unittest.TestCase):

    def testPoisson(self):

        poisson = """((Jinv_00*FE0_D10_ip_j + Jinv_10*FE0_D01_ip_j)*(Jinv_00*FE0_D10_ip_k + Jinv_10*FE0_D01_ip_k) + (Jinv_01*FE0_D10_ip_j + Jinv_11*FE0_D01_ip_j)*(Jinv_01*FE0_D10_ip_k + Jinv_11*FE0_D01_ip_k))*W4_ip*det"""

        expr = Product([
                     Sum([
                          Product([
                                   Sum([
                                        Product([Symbol("Jinv_00", GEO), Symbol("FE0_D10_ip_j", BASIS)])
                                        ,
                                        Product([Symbol("Jinv_10", GEO), Symbol("FE0_D01_ip_j", BASIS)])
                                       ]),
                                   Sum([
                                        Product([Symbol("Jinv_00", GEO), Symbol("FE0_D10_ip_k", BASIS)])
                                        ,
                                        Product([Symbol("Jinv_10", GEO), Symbol("FE0_D01_ip_k", BASIS)])
                                       ])
                                  ])
                          ,
                          Product([
                                   Sum([
                                        Product([Symbol("Jinv_01", GEO), Symbol("FE0_D10_ip_j", BASIS)])
                                        ,
                                        Product([Symbol("Jinv_11", GEO), Symbol("FE0_D01_ip_j", BASIS)])
                                       ]),
                                   Sum([
                                        Product([Symbol("Jinv_01", GEO), Symbol("FE0_D10_ip_k", BASIS)])
                                        ,
                                        Product([Symbol("Jinv_11", GEO), Symbol("FE0_D01_ip_k", BASIS)])
                                       ])
                                  ])
                         ])
                     ,
                     Symbol("W4_ip", IP)
                     ,
                     Symbol("det", GEO)
                    ])

#        print "\nPoisson"
#        start = time.time()
        expr_exp = expr.expand()
#        print "Poisson: time, expand():     ", time.time() - start

#        start = time.time()
        poisson_exp = expand_operations(poisson, get_format())
#        print "Poisson: time, old expand(): ", time.time() - start

#        start = time.time()
        expr_red = expr_exp.reduce_ops()
#        print "Poisson: time, reduce_ops(): ", time.time() - start

#        start = time.time()
        poisson_red = reduce_operations(poisson, get_format())
#        print "Poisson: time, old reduce(): ", time.time() - start

        poisson_exp_ops = operation_count(poisson_exp, get_format())
        poisson_red_ops = operation_count(poisson_red, get_format())
#        print "expr.ops():           ", expr.ops()
#        print "Poisson old exp: ops: ", poisson_exp_ops
#        print "expr_exp.ops():       ", expr_exp.ops()
#        print "Poisson old red: ops: ", poisson_red_ops
#        print "expr_red.ops():       ", expr_red.ops()

#        print "expr: ", expr
#        print "exp:  ", expr_exp
#        print "red:  ", expr_red

        Jinv_00, Jinv_01, Jinv_10, Jinv_11, W4_ip, det = (1.1, 1.5, -4.3, 1.7, 11, 52.3)
        FE0_D01_ip_j, FE0_D10_ip_j, FE0_D01_ip_k, FE0_D10_ip_k = (1.12, 5.7, -9.3, 7.4)
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(poisson)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(poisson_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(poisson_red)))
        self.assertEqual(expr.ops(), 17)
        self.assertEqual(poisson_exp_ops, 47)
        self.assertEqual(expr_exp.ops(), 47)
        self.assertEqual(poisson_red_ops, 23)
        self.assertEqual(expr_red.ops(), 23)

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestPoisson('testPoisson'))

