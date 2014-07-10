#!/usr/bin/env python
"Test suite for the symbolic classes."

# Copyright (C) 2009-2010 Kristian B. Oelgaard
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
# First added:  2009-07-11
# Last changed: 2010-02-01

import unittest

# FFC modules
from ffc.quadrature.symbolics import *

# Import tests
from .testfloat import TestFloat
from .testsymbol import TestSymbol
from .testproduct import TestProduct
from .testsum import TestSum
from .testfraction import TestFraction
from .testfloatoperators import TestFloatOperators
from .testsymboloperators import TestSymbolOperators
from .testproductoperators import TestProductOperators
from .testsumoperators import TestSumOperators
from .testfractionoperators import TestFractionOperators
from .testmixedsymbols import TestMixedSymbols
from .testexpandoperations import TestExpandOperations
from .testreducevartype import TestReduceVarType
from .testreduceoperations import TestReduceOperations
from .testnotfinished import TestNotFinished
from .testdgelastodyn import TestDGElastoDyn
from .testreducegip import TestReduceGIP
from .testpoisson import TestPoisson
from .testelasticity2d import TestElasticity2D
from .testelasticityterm import TestElasticityTerm
from .testelasweighted import TestElasWeighted
from .testelasweighted2 import TestElasWeighted2
from .testrealexamples import TestRealExamples

class TestSingle(unittest.TestCase):

    def testSingle(self):
        "Run a single test."
        expr =\
Fraction(
  Symbol('W1', GEO),
  Sum([
    Fraction(
      Sum([
        Symbol('F1', IP),
        Symbol('F2', IP)
      ]),
      Sum([
        Symbol('K_00', GEO), Symbol('K_01', GEO)
      ])
    ),
    Fraction(
      Sum([
        Symbol('F3', IP),
        Symbol('F4', IP)
      ]),
      Sum([
        Symbol('K_10', GEO), Symbol('K_11', GEO),
      ])
    )
  ])
)
#        print "RED: ", expr
        red = expr.expand().reduce_vartype(IP)
        red = expr.reduce_vartype(IP)
#        red = expr.reduce_vartype(IP)
#        print "expr: ", expr.expand()
#        comb = Sum([Product([f, r]) for f,r in red]).expand()
#        print "\nsum: \n", comb
#        print "eval expr: ", eval(str(expr))
#        print "eval comb: ", eval(str(comb))
#        self.assertAlmostEqual(eval(str(expr)), eval(str(comb)))

def suite():

    suite = unittest.TestSuite()
    # Classes and member functions
    suite.addTest(TestFloat('testFloat'))
    suite.addTest(TestSymbol('testSymbol'))
    suite.addTest(TestProduct('testProduct'))
    suite.addTest(TestSum('testSum'))
    suite.addTest(TestFraction('testFraction'))
    suite.addTest(TestFloatOperators('testFloatOperators'))
    suite.addTest(TestSymbolOperators('testSymbolOperators'))
    suite.addTest(TestProductOperators('testProductOperators'))
    suite.addTest(TestSumOperators('testSumOperators'))
    suite.addTest(TestFractionOperators('testFractionOperators'))
    suite.addTest(TestMixedSymbols('testMixedSymbols'))
    suite.addTest(TestExpandOperations('testExpandOperations'))
    suite.addTest(TestReduceVarType('testReduceVarType'))
    suite.addTest(TestReduceOperations('testReduceOperations'))

    # Misc.
    suite.addTest(TestNotFinished('testNotFinished'))

    # 'Real' expressions (expand and reduce)
    suite.addTest(TestDGElastoDyn('testDGElastoDyn'))
    suite.addTest(TestReduceGIP('testReduceGIP'))
    suite.addTest(TestPoisson('testPoisson'))
    suite.addTest(TestElasticity2D('testElasticity2D'))

    # 'Real' expressions (generate code)
    suite.addTest(TestElasticityTerm('testElasticityTerm'))
    suite.addTest(TestElasWeighted('testElasWeighted'))
    suite.addTest(TestElasWeighted2('testElasWeighted2'))

    # Various bug encounters
    suite.addTest(TestRealExamples('testRealExamples'))

    return suite

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(suite())
#    runner.run(TestSingle('testSingle'))

