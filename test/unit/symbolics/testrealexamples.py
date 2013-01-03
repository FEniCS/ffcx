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

class TestRealExamples(unittest.TestCase):

    def testRealExamples(self):

#        p = Product([
#                    Sum([
#                        Product([
#                                  Symbol('w[5][0]', GEO),
#                                  Fraction(
#                                            Product([
#                                                    Symbol('FE0_C1_D01[ip][k]', BASIS), Symbol('Jinv_10', GEO)
#                                                    ]),
#                                            Product([
#                                                    Symbol('w[5][0]', GEO), Symbol('w[5][0]', GEO)
#                                                    ])
#                                            )
#                                                    
#                                ]),
#                        Product([
#                                  Symbol('w[5][0]', GEO),
#                                  Fraction(
#                                          Product([
#                                                    Symbol('FE0_C1_D01[ip][k]', BASIS), Symbol('Jinv_11', GEO)
#                                                  ]),
#                                          Product([
#                                                    Symbol('w[5][0]', GEO), Symbol('w[5][0]', GEO)
#                                                  ])
#                                          )
#                                ])
#                        ])
#                   ])

#        p = Product([
#                      Sum([
#                            Product([
#                                      Symbol('x', BASIS),
#                                      Sum([
#                                            Symbol('y', BASIS),
#                                            Product([
#                                                      Sum([
#                                                           Symbol('y', BASIS),
#                                                           Product([
#                                                                    Symbol('y', BASIS),
#                                                                    Symbol('z', GEO)
#                                                                   ]),
#                                                           Symbol('y', BASIS)
#                                                         ])
#                                                    ]),
#                                           Symbol('y', BASIS)
#                                          ])
#                                    ]),
#                          Product([
#                                  Symbol('x', BASIS),
#                                  Sum([
#                                        Product([
#                                                Symbol('y', BASIS),
#                                                              Symbol('z', GEO)
#                                              ]),
#                                        Symbol('y', BASIS)
#                                      ])
#                                ])

#                          ])
#                      ])

#        p = Product([
#                     Sum([
#                          Product([
#                                    Symbol('FE0_C1_D01[ip][j]', BASIS),
#                                    Product([
#                                            Symbol('FE0_C1_D01[ip][k]', BASIS),
#                                            Sum([
#                                                 Symbol('w[4][0]', GEO)
#                                                ]),
#                                            Sum([
#                                                  Symbol('w[4][0]', GEO)
#                                                ])
#                                          ])
#                                    ]),
#                        Product([
#                                  Symbol('FE0_C1_D01[ip][j]', BASIS),
#                                  Symbol('FE0_C1_D01[ip][k]', BASIS)
#                                ])
#                         ])
#                      ])

        p = Product([ Symbol('FE0_C1_D01[ip][k]', BASIS),
                      Sum([
                            Symbol('Jinv_10', GEO),
                            Symbol('w[4][0]', GEO)
                          ]),
                      Sum([
                            Symbol('Jinv_10', GEO),
                            Symbol('w[4][0]', GEO)
                          ])
                    ])

#        print "p: ", p
#        print p.expand()

        br = p.reduce_vartype(BASIS)
#        print
#        print br[0]
#        print br[1]

        be = p.expand().reduce_vartype(BASIS)
#        print
#        print be[0][0]
#        print be[0][1]
        if len(be) == 1:
            if be[0][0] == br[0]:
                if be[0][1] != br[1].expand():
#                        print "\np: ", repr(p)
                        print("\nbe: ", repr(be[0][1]))
                        print("\nbr: ", repr(br[1].expand()))
                        print("\nbe: ", be[0][1])
                        print("\nbr: ", br[1].expand())
                        error("here1")

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestRealExamples('testRealExamples'))

