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
from ffc import default_parameters
set_float_formatting(default_parameters()['precision'])

class TestElasWeighted(unittest.TestCase):

    def testElasWeighted(self):
        expr = Product([
                          Symbol('W4', IP),
                          Sum([
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', BASIS),
                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                        Symbol('Jinv_00', GEO),
                                        Symbol('w1', GEO)
                                        ]),
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', BASIS),
                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                        Symbol('Jinv_01', GEO),
                                        Symbol('w0', GEO)
                                        ]),
                              Product([
                                        Symbol('w2', GEO),
                                        Sum([
                                              Product([
                                                      Symbol('FE0_C1_D01_ip_j', BASIS),
                                                      Symbol('FE0_C1_D01_ip_k', BASIS),
                                                      Symbol('Jinv_00', GEO),
                                                      Symbol('w1', GEO)
                                                      ]),
                                              Product([
                                                      Symbol('FE0_C1_D01_ip_j', BASIS),
                                                      Symbol('FE0_C1_D01_ip_k', BASIS),
                                                      Symbol('Jinv_01', GEO),
                                                      Symbol('w0', GEO)
                                                      ])
                                            ])
                                      ])
                              ])
                          ])
                                                       
#        print "\nElasticityWeighted"
#        start = time.time()
        expr_exp = expr.expand()
#        print "ElasWeighted: time, expand():     ", time.time() - start

#        start = time.time()
        expr_red = expr_exp.reduce_ops()
#        print "ElasWeighted: time, reduce_ops(): ", time.time() - start

#        print "expr.ops():     ", expr.ops()
#        print "expr_exp.ops(): ", expr_exp.ops()
#        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        det, W4, w0, w1, w2, Jinv_00, Jinv_01, Jinv_11, Jinv_10, FE0_C1_D01_ip_j, FE0_C1_D01_ip_k = [0.123 + i for i in range(11)]
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertEqual(expr.ops(), 17)
        self.assertEqual(expr_exp.ops(), 21)
        self.assertEqual(expr_red.ops(), 10)

        # Generate code
        ip_consts = {}
        geo_consts = {}
        trans_set = set()

        start = time.time()
        opt_code = optimise_code(expr, ip_consts, geo_consts, trans_set)
#        print "ElasWeighted, optimise_code(): ", time.time() - start

        G0 = eval(str(geo_consts.items()[0][0]))
        Gip0 = eval(str(ip_consts.items()[0][0]))
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestElasWeighted('testElasWeighted'))

