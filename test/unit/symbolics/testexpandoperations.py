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

class TestExpandOperations(unittest.TestCase):

    def testExpandOperations(self):
        f0 = FloatValue(-1)
        f1 = FloatValue(2)
        f2 = FloatValue(1)
        sx = Symbol("x", GEO)
        sy = Symbol("y", GEO)
        sz = Symbol("z", GEO)
        s0 = Product([FloatValue(-1), Symbol("x", GEO)])
        s1 = Symbol("y", GEO)
        s2 = Product([FloatValue(5), Symbol("z", IP)])
        s3 = Product([FloatValue(-4), Symbol("z", GEO)])

        # Random variable values
        x = 2.2
        y = -0.2
        z = 1.1

        # Aux. expressions
        P0 = Product([s2, s1])
        P1 = Product([P0, s0])
        P2 = Product([P1, s1, P0])
        P3 = Product([P1, P2])

        S0 = Sum([s2, s1])
        S1 = Sum([S0, s0])
        S2 = Sum([S1, s1, S0])
        S3 = Sum([S1, S2])

        F0 = Fraction(s2, s1)
        F1 = Fraction(F0, s0)
        F2 = Fraction(F1, F0)
        F3 = Fraction(F1, F2)

        # Special fractions
        F4 = Fraction(P0, F0)
        F5 = Fraction(Fraction(s0, P0), P0)
        F6 = Fraction( Fraction( Fraction(s1, s0), Fraction(s1, s2)), Fraction( Fraction(s2, s0), Fraction(s1, s0)) )
        F7 = Fraction(s1, Product([s1, Symbol("x", GEO)]))
        F8 = Fraction( Sum([sx, Fraction(sy, sx)]), FloatValue(2))

        F4x = F4.expand()
        F5x = F5.expand()
        F6x = F6.expand()
        F7x = F7.expand()
        F8x = F8.expand()

#        print "\nF4: '%s'" %F4
#        print "F4x: '%s'" %F4x
#        print "\nF5: '%s'" %F5
#        print "F5x: '%s'" %F5x
#        print "\nF6: '%s'" %F6
#        print "F6x: '%s'" %F6x
#        print "\nF7: '%s'" %F7
#        print "F7x: '%s'" %F7x
#        print "\nF8: '%s'" %F8
#        print "F8x: '%s'" %F8x

        self.assertAlmostEqual(eval(str(F4)), eval(str(F4x)))
        self.assertAlmostEqual(eval(str(F5)), eval(str(F5x)))
        self.assertAlmostEqual(eval(str(F6)), eval(str(F6x)))
        self.assertAlmostEqual(eval(str(F7)), eval(str(F7x)))
        self.assertAlmostEqual(eval(str(F8)), eval(str(F8x)))

        self.assertEqual(F4.ops(), 5)
        self.assertEqual(F4x.ops(), 1)
        self.assertEqual(F5.ops(), 6)
        self.assertEqual(F5x.ops(), 5)
        self.assertEqual(F6.ops(), 9)
        self.assertEqual(F6x.ops(), 1)
        self.assertEqual(F7.ops(), 2)
        self.assertEqual(F7x.ops(), 1)
        self.assertEqual(F8.ops(), 3)
        self.assertEqual(F8x.ops(), 4)

        # Expressions that should be expanded
        e0 = Product([P3, F2])
        e1 = Product([S3, P2])
        e2 = Product([F3, S1])

        e3 = Sum([P3, F2])
        e4 = Sum([S3, P2])
        e5 = Sum([F3, S1])

        e6 = Fraction(P3, F2)
        e7 = Fraction(S3, P2)
        e8 = Fraction(F3, S1)
        e9 = Fraction(S0, s0)

        e0x = e0.expand()
        e1x = e1.expand()
        e2x = e2.expand()
        e3x = e3.expand()
        e4x = e4.expand()
        e5x = e5.expand()
        e6x = e6.expand()
        e7x = e7.expand()
        e8x = e8.expand()
        e9x = e9.expand()

#        print "\ne0: '%s'" %e0
#        print "e0x: '%s'" %e0x
#        print "\ne1: '%s'" %e1
#        print "e1x: '%s'" %e1x
#        print "\ne2: '%s'" %e2
#        print "e2x: '%s'" %e2x
#        print "\ne3: '%s'" %e3
#        print "e3x: '%s'" %e3x
#        print "\ne4: '%s'" %e4
#        print "e4x: '%s'" %e4x
#        print "\ne5: '%s'" %e5
#        print "e5x: '%s'" %e5x
#        print "\ne6: '%s'" %e6
#        print "e6x: '%s'" %e6x
#        print "\ne7: '%s'" %e7
#        print "e7x: '%s'" %e7x
#        print "\ne8: '%s'" %e8
#        print "e8x: '%s'" %e8x
#        print "\ne9: '%s'" %e9
#        print "e9x: '%s'" %e9x

        self.assertAlmostEqual(eval(str(e0)), eval(str(e0x)))
        self.assertAlmostEqual(eval(str(e1)), eval(str(e1x)))
        self.assertAlmostEqual(eval(str(e2)), eval(str(e2x)))
        self.assertAlmostEqual(eval(str(e3)), eval(str(e3x)))
        self.assertAlmostEqual(eval(str(e4)), eval(str(e4x)))
        self.assertAlmostEqual(eval(str(e5)), eval(str(e5x)))
        self.assertAlmostEqual(eval(str(e6)), eval(str(e6x)))
        self.assertAlmostEqual(eval(str(e7)), eval(str(e7x)))
        self.assertAlmostEqual(eval(str(e8)), eval(str(e8x)))
        self.assertAlmostEqual(eval(str(e9)), eval(str(e9x)))

        self.assertEqual(e0.ops(), 16)
        self.assertEqual(e0x.ops(), 8)
        self.assertEqual(e1.ops(), 18)
        self.assertEqual(e1x.ops(), 23)
        self.assertEqual(e2.ops(), 14)
        self.assertEqual(e2x.ops(), 9)
        self.assertEqual(e3.ops(), 16)
        self.assertEqual(e3x.ops(), 11)
        self.assertEqual(e4.ops(), 18)
        self.assertEqual(e4x.ops(), 12)
        self.assertEqual(e5.ops(), 14)
        self.assertEqual(e5x.ops(), 6)
        self.assertEqual(e6.ops(), 16)
        self.assertEqual(e6x.ops(), 10)
        self.assertEqual(e7.ops(), 18)
        self.assertEqual(e7x.ops(), 17)
        self.assertEqual(e8.ops(), 14)
        self.assertEqual(e8x.ops(), 8)
        self.assertEqual(e9.ops(), 3)
        self.assertEqual(e9x.ops(), 4)

        # More expressions (from old expand tests)
        PF = Product([F0, F1])
        E0 = Product([s1, f0, S1])
        E1 = Sum([P0, E0])
        E2 = Fraction(Sum([Product([f1])]), f2)
        E3 = Sum([F0, F0])
        E4 = Product([ Sum([ Product([sx, Sum([sy, Product([ Sum([sy, Product([sy, sz]), sy])]), sy])]),
                             Product([sx, Sum([ Product([sy, sz]), sy])])])])
        P4 = Product([s1,
        Sum([s0, s1])])
        P5 = Product([s0, E0])
        P6 = Product([s1])
        S4 = Sum([s1])


        # Create 'real' term that caused me trouble
        P00 = Product([Symbol("Jinv_00", GEO)]*2)
        P01 = Product([Symbol("Jinv_01", GEO)]*2)
        P20 = Product([Symbol("Jinv_00", GEO),
        Product([f1, Symbol("Jinv_20", GEO)]) ])
        P21 = Product([Symbol("Jinv_01", GEO),
        Product([f1, Symbol("Jinv_21", GEO)]) ])
        PS0 = Product([Symbol("Jinv_22", GEO),
        Sum([P00, P01])])
        PS1 = Product([ Product([f0, Symbol("Jinv_02", GEO)]),
        Sum([P20, P21])])
        SP = Sum([PS0, PS1])

        PFx = PF.expand()
        E0x = E0.expand()
        E1x = E1.expand()
        E2x = E2.expand()
        E3x = E3.expand()
        E4x = E4.expand()
        P4x = P4.expand()
        P5x = P5.expand()
        P6x = P6.expand()
        S4x = S4.expand()
        SPx = SP.expand()

#        print "\nPF: '%s'" %PF
#        print "PFx: '%s'" %PFx
#        print "\nE0: '%s'" %E0
#        print "E0x: '%s'" %E0x
#        print "\nE1: '%s'" %E1
#        print "E1x: '%s'" %E1x
#        print "\nE2: '%s'" %E2
#        print "E2x: '%s'" %E2x
#        print "\nE3: '%s'" %E3
#        print "E3x: '%s'" %E3x
#        print "\nE4: '%s'" %E4
#        print "E4x: '%s'" %E4x
#        print "\nP4: '%s'" %P4
#        print "P4x: '%s'" %P4x
#        print "\nP5: '%s'" %P5
#        print "P5x: '%s'" %P5x
#        print "\nP6: '%s'" %repr(P6)
#        print "P6x: '%s'" %repr(P6x)
#        print "\nS4: '%s'" %repr(S4)
#        print "S4x: '%s'" %repr(S4x)
#        print "\nSP: '%s'" %SP
#        print "SPx: '%s'" %SPx

        Jinv_00, Jinv_01, Jinv_10, Jinv_02, Jinv_20, Jinv_22, Jinv_21, W1, det = [1,2,3,4,5,6,7,8,9]

        self.assertAlmostEqual(eval(str(SP)), eval(str(SPx)))
        self.assertAlmostEqual(eval(str(E0)), eval(str(E0x)))
        self.assertAlmostEqual(eval(str(E1)), eval(str(E1x)))
        self.assertAlmostEqual(eval(str(E2)), eval(str(E2x)))
        self.assertAlmostEqual(eval(str(E3)), eval(str(E3x)))
        self.assertAlmostEqual(eval(str(E4)), eval(str(E4x)))
        self.assertAlmostEqual(eval(str(SP)), eval(str(SPx)))
        self.assertAlmostEqual(eval(str(P4)), eval(str(P4x)))
        self.assertAlmostEqual(eval(str(P5)), eval(str(P5x)))
        self.assertEqual(P6x, s1)
        self.assertEqual(S4x, s1)
        self.assertEqual(PF.ops(), 6)
        self.assertEqual(PFx.ops(), 5)
        self.assertEqual(E0.ops(), 4)
        self.assertEqual(E0x.ops(), 6)
        self.assertEqual(E1.ops(), 7)
        self.assertEqual(E1x.ops(), 3)
        self.assertEqual(E2.ops(), 1)
        self.assertEqual(E2x.ops(), 0)
        self.assertEqual(E3.ops(), 5)
        self.assertEqual(E3x.ops(), 5)
        self.assertEqual(E4.ops(), 10)
        self.assertEqual(E4x.ops(), 6)
        self.assertEqual(SP.ops(), 11)
        self.assertEqual(SPx.ops(), 13)
        self.assertEqual(P4.ops(), 2)
        self.assertEqual(P4x.ops(), 3)
        self.assertEqual(P5.ops(), 5)
        self.assertEqual(P5x.ops(), 9)

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestExpandOperations('testExpandOperations'))

