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
# Last changed: 2010-03-11

# Pyhton modules
import unittest
import time

# FFC modules
from ffc.quadrature.symbolics import *
from ffc.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])

class TestMixedSymbols(unittest.TestCase):

    def testMixedSymbols(self):

        f_0 = format["float"](0)
        f_2 = format["float"](2)
        f_3 = format["float"](3)
        f_4 = format["float"](4)
        f_6 = format["float"](6)

        f0 = FloatValue(-2.0)
        f1 = FloatValue(3.0)
        f2 = FloatValue(0)

        s0 = Symbol("x", BASIS)
        s1 = Symbol("y", GEO)
        s2 = Symbol("z", GEO)

        p0 = Product([s0, s1])
        p1 = Product([f1, s0, s1])
        p2 = Product([s0, f2, s2])
        p3 = Product([s0, f0, s1, f1, s2])

        S0 = Sum([s0, s1])
        S1 = Sum([s0, s0])
        S2 = Sum([f0, s0])
        S3 = Sum([s0, f0, s0])

        F0 = Fraction(f1, f0)
        F1 = Fraction(s0, s1)
        F2 = Fraction(s0, f1)
        F3 = Fraction(f0, s1)

        x = 1.2; y = 2.36; z = 6.75;
        # Mixed products
        mpp0 = Product([p0, s0])
        mpp1 = Product([p1, p0])
        mpp2 = Product([p2, p3])
        mpp3 = Product([p1, mpp1])

        mps0 = Product([S0, s0])
        mps1 = Product([S1, S0])
        mps2 = Product([S2, S3])
        mps3 = Product([S1, mps1])

        mpf0 = Product([F1, s0])
        mpf1 = Product([F1, F2])
        mpf2 = Product([F2, F3])
        mpf3 = Product([F1, mpf1])

#        print "\nMixed Products"
#        print "\nmpp0: %s * %s = '%s'" % (p0, s0, mpp0)
#        print "mpp1: %s * %s = '%s'" % (p1, p0, mpp1)
#        print "mpp2: %s * %s = '%s'" % (p2, p3, mpp2)
#        print "mpp3: %s * %s = '%s'" % (p1, mpp1, mpp3)
#        print "\nmps0: %s * %s = '%s'" % (S0, s0, mps0)
#        print "mps1: %s * %s = '%s'" % (S1, S0, mps1)
#        print "mps2: %s * %s = '%s'" % (S2, S3, mps2)
#        print "mps3: %s * %s = '%s'" % (S1, mps1, mps3)
#        print "\nmpf0: %s * %s = '%s'" % (F1, s0, mpf0)
#        print "mpf1: %s * %s = '%s'" % (F1, F2, mpf1)
#        print "mpf2: %s * %s = '%s'" % (F2, F3, mpf2)
#        print "mpf3: %s * %s = '%s'" % (F1, mpf1, mpf3)

        self.assertAlmostEqual(eval(str(mpp0)), eval(str(p0))*eval(str(s0)))
        self.assertAlmostEqual(eval(str(mpp1)), eval(str(p1))*eval(str(p0)))
        self.assertAlmostEqual(eval(str(mpp2)), eval(str(p2))*eval(str(p3)))
        self.assertAlmostEqual(eval(str(mpp3)), eval(str(p1))*eval(str(mpp1)))

        self.assertAlmostEqual(eval(str(mps0)), eval(str(S0))*eval(str(s0)))
        self.assertAlmostEqual(eval(str(mps1)), eval(str(S1))*eval(str(S0)))
        self.assertAlmostEqual(eval(str(mps2)), eval(str(S2))*eval(str(S3)))
        self.assertAlmostEqual(eval(str(mps3)), eval(str(S1))*eval(str(mps1)))

        self.assertAlmostEqual(eval(str(mpf0)), eval(str(F1))*eval(str(s0)))
        self.assertAlmostEqual(eval(str(mpf1)), eval(str(F1))*eval(str(F2)))
        self.assertAlmostEqual(eval(str(mpf2)), eval(str(F2))*eval(str(F3)))
        self.assertAlmostEqual(eval(str(mpf3)), eval(str(F1))*eval(str(mpf1)))

        self.assertEqual(mpp0.ops(), 2)
        self.assertEqual(mpp1.ops(), 4)
        self.assertEqual(mpp2.ops(), 0)
        self.assertEqual(mpp3.ops(), 6)

        self.assertEqual(mps0.ops(), 2)
        self.assertEqual(mps1.ops(), 3)
        self.assertEqual(mps2.ops(), 4)
        self.assertEqual(mps3.ops(), 5)

        self.assertEqual(mpf0.ops(), 2)
        self.assertEqual(mpf1.ops(), 3)
        self.assertEqual(mpf2.ops(), 3)
        self.assertEqual(mpf3.ops(), 5)

        self.assertEqual(str(mpp0), 'x*x*y')
        self.assertEqual(str(mpp1), '%s*x*x*y*y' % f_3)
        self.assertEqual(str(mpp2), '%s' % f_0)
        self.assertEqual(str(mpp3), '%s*x*x*x*y*y*y' % format["float"](9))
        self.assertEqual(str(mps0), 'x*(x + y)')
        self.assertEqual(str(mps1), '(x + x)*(x + y)')
#        self.assertEqual(str(mps2), '(x-2)*(x + x-2)')
        self.assertEqual(str(mps2), '(x + x-%s)*(x-%s)' % (f_2, f_2))
        self.assertEqual(str(mps3), '(x + x)*(x + x)*(x + y)')
        self.assertEqual(str(mpf0), 'x*x/y')
        self.assertEqual(str(mpf1), 'x/%s*x/y' % f_3)
        self.assertEqual(str(mpf2), '-%s/y*x/%s' % (f_2, f_3))
        self.assertEqual(str(mpf3), 'x/%s*x/y*x/y' % f_3)


        # Mixed sums
        msp0 = Sum([p0, s0])
        msp1 = Sum([p1, p0])
        msp2 = Sum([p2, p3])
        msp3 = Sum([p1, msp1])
        msp4 = Sum([f2, f2])

        mss0 = Sum([S0, s0])
        mss1 = Sum([S1, S0])
        mss2 = Sum([S2, S3])
        mss3 = Sum([S1, mps1])

        msf0 = Sum([F1, s0])
        msf1 = Sum([F1, F2])
        msf2 = Sum([F2, F3])
        msf3 = Sum([F1, msf1])

#        print "\nTesting Mixed Sums"
#        print "\nmsp0: %s + %s = '%s'" % (p0, s0, msp0)
#        print "msp1: %s + %s = '%s'" % (p1, p0, msp1)
#        print "msp2: %s + %s = '%s'" % (p2, p3, msp2)
#        print "msp3: %s + %s = '%s'" % (p1, msp1, msp3)
#        print "msp4: %s + %s = '%s'" % (f2, f2, msp4)
#        print "\nmss0: %s + %s = '%s'" % (S0, s0, mss0)
#        print "mss1: %s + %s = '%s'" % (S1, S0, mss1)
#        print "mss2: %s + %s = '%s'" % (S2, S3, mss2)
#        print "mss3: %s + %s = '%s'" % (S1, mss1, mss3)
#        print "\nmsf0: %s + %s = '%s'" % (F1, s0, msf0)
#        print "msf1: %s + %s = '%s'" % (F1, F2, msf1)
#        print "msf2: %s + %s = '%s'" % (F2, F3, msf2)
#        print "msf3: %s + %s = '%s'" % (F1, msf1, msf3)
#        print "msf3: %s + %s = '%s'" % (F1, msf1, msf3)

        self.assertAlmostEqual(eval(str(msp0)), eval(str(p0))+eval(str(s0)))
        self.assertAlmostEqual(eval(str(msp1)), eval(str(p1))+eval(str(p0)))
        self.assertAlmostEqual(eval(str(msp2)), eval(str(p2))+eval(str(p3)))
        self.assertAlmostEqual(eval(str(msp3)), eval(str(p1))+eval(str(msp1)))
        self.assertEqual(str(msp4), '%s' % f_0)

        self.assertAlmostEqual(eval(str(mss0)), eval(str(S0))+eval(str(s0)))
        self.assertAlmostEqual(eval(str(mss1)), eval(str(S1))+eval(str(S0)))
        self.assertAlmostEqual(eval(str(mss2)), eval(str(S2))+eval(str(S3)))
        self.assertAlmostEqual(eval(str(mss3)), eval(str(S1))+eval(str(mps1)))

        self.assertAlmostEqual(eval(str(msf0)), eval(str(F1))+eval(str(s0)))
        self.assertAlmostEqual(eval(str(msf1)), eval(str(F1))+eval(str(F2)))
        self.assertAlmostEqual(eval(str(msf2)), eval(str(F2))+eval(str(F3)))
        self.assertAlmostEqual(eval(str(msf3)), eval(str(F1))+eval(str(msf1)))

        self.assertEqual(msp0.ops(), 2)
        self.assertEqual(msp1.ops(), 4)
        self.assertEqual(msp2.ops(), 3)
        self.assertEqual(msp3.ops(), 7)

        self.assertEqual(mss0.ops(), 2)
        self.assertEqual(mss1.ops(), 3)
        self.assertEqual(mss2.ops(), 3)
        self.assertEqual(mss3.ops(), 5)

        self.assertEqual(msf0.ops(), 2)
        self.assertEqual(msf1.ops(), 3)
        self.assertEqual(msf2.ops(), 3)
        self.assertEqual(msf3.ops(), 5)

        self.assertEqual(str(msp0), '(x + x*y)')
        self.assertEqual(str(msp1), '(%s*x*y + x*y)' % f_3)
        self.assertEqual(str(msp2), '-%s*x*y*z' % f_6)
        self.assertEqual(str(msp3), '(%s*x*y + %s*x*y + x*y)' % (f_3, f_3))
        self.assertEqual(str(mss0), '(x + x + y)')
        self.assertEqual(str(mss1), '(x + x + x + y)')
        self.assertEqual(str(mss2), '(x + x + x-%s)' % f_4)
        self.assertEqual(str(mss3), '(x + x + (x + x)*(x + y))')
        self.assertEqual(str(msf0), '(x + x/y)')
        self.assertEqual(str(msf1), '(x/%s + x/y)' % f_3)
        self.assertEqual(str(msf2), '(x/%s-%s/y)' % (f_3, f_2))
        self.assertEqual(str(msf3), '(x/%s + x/y + x/y)' % f_3)


        # Mixed fractions
        mfp0 = Fraction(p0, s0)
        mfp1 = Fraction(p1, p0)
        mfp2 = Fraction(p2, p3)
        mfp3 = Fraction(p1, mfp1)

        mfs0 = Fraction(S0, s0)
        mfs1 = Fraction(S1, S0)
        mfs2 = Fraction(S2, S3)
        mfs3 = Fraction(S1, mfs1)

        mff0 = Fraction(F1, s0)
        mff1 = Fraction(F1, F2)
        mff2 = Fraction(F2, F3)
        mff3 = Fraction(F1, mff1)

#        print "\nTesting Mixed Fractions"
#        print "\nmfp0: %s / %s = '%s'" % (p0, s0, mfp0)
#        print "mfp1: %s / %s = '%s'" % (p1, p0, mfp1)
#        print "mfp2: %s / %s = '%s'" % (p2, p3, mfp2)
#        print "mfp3: %s / %s = '%s'" % (p1, mfp1, mfp3)
#        print "\nmfs0: %s / %s = '%s'" % (S0, s0, mfs0)
#        print "mfs1: %s / %s = '%s'" % (S1, S0, mfs1)
#        print "mfs2: %s / %s = '%s'" % (S2, S3, mfs2)
#        print "mfs3: %s / %s = '%s'" % (S1, mfs1, mfs3)
#        print "\nmff0: %s / %s = '%s'" % (F1, s0, mff0)
#        print "mff1: %s / %s = '%s'" % (F1, F2, mff1)
#        print "mff2: %s / %s = '%s'" % (F2, F3, mff2)
#        print "mff3: %s / %s = '%s'" % (F1, mff1, mff3)

        self.assertAlmostEqual(eval(str(mfp0)), eval(str(p0))/eval(str(s0)))
        self.assertAlmostEqual(eval(str(mfp1)), eval(str(p1))/eval(str(p0)))
        self.assertAlmostEqual(eval(str(mfp2)), eval(str(p2))/eval(str(p3)))
        self.assertAlmostEqual(eval(str(mfp3)), eval(str(p1))/eval(str(mfp1)))

        self.assertAlmostEqual(eval(str(mfs0)), eval(str(S0))/eval(str(s0)))
        self.assertAlmostEqual(eval(str(mfs1)), eval(str(S1))/eval(str(S0)))
        self.assertAlmostEqual(eval(str(mfs2)), eval(str(S2))/eval(str(S3)))
        self.assertAlmostEqual(eval(str(mfs3)), eval(str(S1))/eval(str(mfs1)))

        self.assertAlmostEqual(eval(str(mff0)), eval(str(F1))/eval(str(s0)))
        self.assertAlmostEqual(eval(str(mff1)), eval(str(F1))/eval(str(F2)))
        self.assertAlmostEqual(eval(str(mff2)), eval(str(F2))/eval(str(F3)))
        self.assertAlmostEqual(eval(str(mff3)), eval(str(F1))/eval(str(mff1)))

        self.assertEqual(mfp0.ops(), 2)
        self.assertEqual(mfp1.ops(), 4)
        self.assertEqual(mfp2.ops(), 0)
        self.assertEqual(mfp3.ops(), 7)

        self.assertEqual(mfs0.ops(), 2)
        self.assertEqual(mfs1.ops(), 3)
        self.assertEqual(mfs2.ops(), 4)
        self.assertEqual(mfs3.ops(), 5)

        self.assertEqual(mff0.ops(), 2)
        self.assertEqual(mff1.ops(), 3)
        self.assertEqual(mff2.ops(), 3)
        self.assertEqual(mff3.ops(), 5)

        self.assertEqual(str(mfp0), 'x*y/x')
        self.assertEqual(str(mfp1), '%s*x*y/(x*y)' % f_3)
        self.assertEqual(str(mfp2), '%s' % f_0)
        self.assertEqual(str(mfp3), '%s*x*y/(%s*x*y/(x*y))' % (f_3, f_3))
        self.assertEqual(str(mfs0), '(x + y)/x')
        self.assertEqual(str(mfs1), '(x + x)/(x + y)')
        self.assertEqual(str(mfs2), '(x-%s)/(x + x-%s)' % (f_2, f_2))
        self.assertEqual(str(mfs3), '(x + x)/((x + x)/(x + y))')
        self.assertEqual(str(mff0), '(x/y)/x')
        self.assertEqual(str(mff1), '(x/y)/(x/%s)' % f_3)
        self.assertEqual(str(mff2), '(x/%s)/(-%s/y)' % (f_3, f_2))
        self.assertEqual(str(mff3), '(x/y)/((x/y)/(x/%s))' % f_3)

        # Use p1 as a base expression for Symbol
        s3 = Symbol(format["cos"], CONST, p1, 1)
        self.assertEqual(str(s3), 'std::cos(%s*x*y)' % f_3)
        self.assertEqual(s3.ops(), 3)

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestMixedSymbols('testMixedSymbols'))

