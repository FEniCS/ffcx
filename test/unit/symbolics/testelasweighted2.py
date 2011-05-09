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
# Last changed: 2010-03-11

# Pyhton modules
import unittest
import time

# FFC modules
from ffc.quadrature.symbolics import *
from ffc.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])

class TestElasWeighted2(unittest.TestCase):

    def testElasWeighted2(self):

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
                                        Symbol('Jinv_01', GEO),
                                        Sum([
                                              Product([
                                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                                        Symbol('w0', GEO)
                                                        ]),
                                              Product([
                                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                                        Symbol('w1', GEO)
                                                        ])
                                              ])
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
                                                    Symbol('Jinv_01', GEO),
                                                    Sum([
                                                          Product([
                                                                  Symbol('FE0_C1_D01_ip_k', BASIS),
                                                                  Symbol('w0', GEO)
                                                                  ]),
                                                          Product([
                                                                   Symbol('FE0_C1_D01_ip_k', BASIS),
                                                                   Symbol('w1', GEO)
                                                                  ])
                                                          ])
                                                    ])
                                            ])
                                        ])
                              ])
                        ])

#        print "\nElasticityWeighted2"
        start = time.time()
        expr_exp = expr.expand()
#        print "ElasWeighted2: time, expand():     ", time.time() - start

#        start = time.time()
        expr_red = expr_exp.reduce_ops()
#        print "ElasWeighted2: time, reduce_ops(): ", time.time() - start

#        print "expr.ops():     ", expr.ops()
#        print "expr_exp.ops(): ", expr_exp.ops()
#        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        det, W4, w0, w1, w2, Jinv_00, Jinv_01, Jinv_11, Jinv_10, FE0_C1_D01_ip_j, FE0_C1_D01_ip_k = [0.123 + i for i in range(11)]
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertEqual(expr.ops(), 21)
        self.assertEqual(expr_exp.ops(), 32)
        self.assertEqual(expr_red.ops(), 13)

        # Generate code
        ip_consts = {}
        geo_consts = {}
        trans_set = set()

        start = time.time()
        opt_code = optimise_code(expr, ip_consts, geo_consts, trans_set)
#        print "ElasWeighted2, optimise_code(): ", time.time() - start

        G = [eval(str(geo_consts.items()[0][0]))]
        I = [eval(str(ip_consts.items()[0][0]))]
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestElasWeighted2('testElasWeighted2'))

