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
from ffc.quadrature.reduce_operations import operation_count, expand_operations, reduce_operations
from ffc.quadrature.symbolics import *
from ffc.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])
class TestElasticity2D(unittest.TestCase):

    def testElasticity2D(self):
        elasticity = """(((Jinv_00*FE0_C0_D10_ip_j + Jinv_10*FE0_C0_D01_ip_j)*2*(Jinv_00*FE0_C0_D10_ip_k + Jinv_10*FE0_C0_D01_ip_k)*2 + ((Jinv_00*FE0_C1_D10_ip_j + Jinv_10*FE0_C1_D01_ip_j) + (Jinv_01*FE0_C0_D10_ip_j + Jinv_11*FE0_C0_D01_ip_j))*((Jinv_00*FE0_C1_D10_ip_k + Jinv_10*FE0_C1_D01_ip_k) + (Jinv_01*FE0_C0_D10_ip_k + Jinv_11*FE0_C0_D01_ip_k))) + ((Jinv_01*FE0_C1_D10_ip_j + Jinv_11*FE0_C1_D01_ip_j)*2*(Jinv_01*FE0_C1_D10_ip_k + Jinv_11*FE0_C1_D01_ip_k)*2 + ((Jinv_01*FE0_C0_D10_ip_j + Jinv_11*FE0_C0_D01_ip_j) + (Jinv_00*FE0_C1_D10_ip_j + Jinv_10*FE0_C1_D01_ip_j))*((Jinv_01*FE0_C0_D10_ip_k + Jinv_11*FE0_C0_D01_ip_k) + (Jinv_00*FE0_C1_D10_ip_k + Jinv_10*FE0_C1_D01_ip_k))))*0.25*W4_ip*det"""

        expr = Product([
                     Sum([
                          Sum([
                               Product([
                                        Sum([
                                             Product([Symbol("Jinv_00", GEO), Symbol("FE0_C0_D10_ip_j", BASIS)])
                                             ,
                                             Product([Symbol("Jinv_10", GEO), Symbol("FE0_C0_D01_ip_j", BASIS)])
                                            ])
                                        ,
                                        FloatValue(2)
                                        ,
                                        Sum([
                                             Product([Symbol("Jinv_00", GEO), Symbol("FE0_C0_D10_ip_k", BASIS)])
                                             ,
                                             Product([Symbol("Jinv_10", GEO), Symbol("FE0_C0_D01_ip_k", BASIS)])
                                            ])
                                        ,
                                        FloatValue(2)
                                        ])
                               ,
                               Product([
                                        Sum([
                                             Sum([
                                                  Product([Symbol("Jinv_00", GEO), Symbol("FE0_C1_D10_ip_j", BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_10", GEO), Symbol("FE0_C1_D01_ip_j", BASIS)])
                                                 ])
                                             ,
                                             Sum([
                                                  Product([Symbol("Jinv_01", GEO), Symbol("FE0_C0_D10_ip_j", BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_11", GEO), Symbol("FE0_C0_D01_ip_j", BASIS)])
                                                 ])
                                            ])
                                        ,
                                        Sum([
                                             Sum([
                                                  Product([Symbol("Jinv_00", GEO), Symbol("FE0_C1_D10_ip_k", BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_10", GEO), Symbol("FE0_C1_D01_ip_k", BASIS)])
                                                 ])
                                             ,
                                             Sum([
                                                  Product([Symbol("Jinv_01", GEO), Symbol("FE0_C0_D10_ip_k", BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_11", GEO), Symbol("FE0_C0_D01_ip_k", BASIS)])
                                                 ])
                                            ])
                                       ])
                              ])
                          ,
                          Sum([
                               Product([
                                        Sum([
                                             Product([Symbol("Jinv_01", GEO), Symbol("FE0_C1_D10_ip_j", BASIS)])
                                             ,
                                             Product([Symbol("Jinv_11", GEO), Symbol("FE0_C1_D01_ip_j", BASIS)])
                                            ])
                                        ,
                                        FloatValue(2)
                                        ,
                                        Sum([
                                             Product([Symbol("Jinv_01", GEO), Symbol("FE0_C1_D10_ip_k", BASIS)])
                                             ,
                                             Product([Symbol("Jinv_11", GEO), Symbol("FE0_C1_D01_ip_k", BASIS)])
                                            ])
                                        ,
                                        FloatValue(2)
                                       ])
                               ,
                               Product([
                                        Sum([
                                             Sum([
                                                  Product([Symbol("Jinv_01", GEO), Symbol("FE0_C0_D10_ip_j", BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_11", GEO), Symbol("FE0_C0_D01_ip_j", BASIS)])
                                                 ])
                                             ,
                                             Sum([
                                                  Product([Symbol("Jinv_00", GEO), Symbol("FE0_C1_D10_ip_j", BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_10", GEO), Symbol("FE0_C1_D01_ip_j", BASIS)])
                                                 ])
                                            ])
                                        ,
                                        Sum([
                                             Sum([
                                                  Product([Symbol("Jinv_01", GEO), Symbol("FE0_C0_D10_ip_k", BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_11", GEO), Symbol("FE0_C0_D01_ip_k", BASIS)])
                                                 ])
                                             ,
                                             Sum([
                                                  Product([Symbol("Jinv_00", GEO), Symbol("FE0_C1_D10_ip_k", BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_10", GEO), Symbol("FE0_C1_D01_ip_k", BASIS)])
                                                 ])
                                            ])
                                       ])
                              ])
                     ])
                     ,
                     FloatValue(0.25)
                     ,
                     Symbol("W4_ip", IP)
                     ,
                     Symbol("det", GEO)
                     ])

#        print "\nElasticity2D"
#        start = time.time()
        expr_exp = expr.expand()
#        print "Elasticity2D: time, expand():     ", time.time() - start

#        start = time.time()
        elasticity_exp = expand_operations(elasticity, format)
#        print "Elasticity2D: time, old expand(): ", time.time() - start

#        start = time.time()
        expr_red = expr_exp.reduce_ops()
#        print "Elasticity2D: time, reduce_ops(): ", time.time() - start

#        start = time.time()
        elasticity_red = reduce_operations(elasticity, format)
#        print "Elasticity2D: time, old reduce(): ", time.time() - start

        elasticity_exp_ops = operation_count(elasticity_exp, format)
        elasticity_red_ops = operation_count(elasticity_red, format)
#        print "expr.ops():                ", expr.ops()
#        print "Elasticity2D old exp: ops: ", elasticity_exp_ops
#        print "expr_exp.ops():            ", expr_exp.ops()
#        print "Elasticity2D old red: ops: ", elasticity_red_ops
#        print "expr_red.ops():            ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red
#        print "old red:\n", elasticity_red

        Jinv_00, Jinv_01, Jinv_10, Jinv_11, W4_ip, det = (1.1, 1.5, -4.3, 1.7, 11, 52.3)
        FE0_C0_D01_ip_j, FE0_C0_D10_ip_j, FE0_C0_D01_ip_k, FE0_C0_D10_ip_k = (1.12, 5.7, -9.3, 7.4)
        FE0_C1_D01_ip_j, FE0_C1_D10_ip_j, FE0_C1_D01_ip_k, FE0_C1_D10_ip_k = (3.12, -8.1, -45.3, 17.5)

        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(elasticity)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(elasticity_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(elasticity_red)))
        self.assertEqual(expr.ops(), 52)
        self.assertEqual(elasticity_exp_ops, 159)
        self.assertEqual(expr_exp.ops(), 159)
        self.assertEqual(elasticity_red_ops, 71)
        self.assertEqual(expr_red.ops(), 71)

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestElasticity2D('testElasticity2D'))

