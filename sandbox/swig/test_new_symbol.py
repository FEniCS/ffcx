#!/usr/bin/env python
"Some simple functions for manIPulating expressions symbolically"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-11 -- 2009-07-15"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"


import unittest
from ufl.common import dstr
from ffc.compiler.codegeneration.quadrature.reduce_operations import operation_count, expand_operations, reduce_operations
from new_symbol import *

import time

from ffc.compiler.format.ufcformat import Format
from ffc.common.constants import FFC_OPTIONS

class Tests(unittest.TestCase):

    def testFloat(self):
        "Test simple FloatValue instance."
        f0 = FloatValue(1.5)
        f1 = FloatValue(-5)
        f2 = FloatValue(-1e-13)
        f3 = FloatValue(-1e-11)
        f4 = FloatValue(1.5)

#        print "\nTesting FloatValue"
#        print "f0: '%s'" %f0
#        print "f1: '%s'" %f1
#        print "f2: '%s'" %f2
#        print "f3: '%s'" %f3

        self.assertEqual(repr(f0), "FloatValue(1.5)")
        self.assertEqual(repr(f1), "FloatValue(-5)")
        self.assertEqual(repr(f2), "FloatValue(0)")
        self.assertEqual(repr(f3), "FloatValue(-1e-11)")

        self.assertEqual(f2.val == 0, True)
        self.assertEqual(f3.val == 0, False)

        self.assertEqual(f0.ops(), 0)
        self.assertEqual(f1.ops(), 0)
        self.assertEqual(f2.ops(), 0)
        self.assertEqual(f3.ops(), 0)

        self.assertEqual(f0 == f4, True)
        self.assertEqual(f1 != f3, True)
        self.assertEqual(f0 < f1, False)
        self.assertEqual(f2 > f3, True)

        # Test hash
        l = [f0]
        d = {f0:0}
        self.assertEqual(f0 in l, True)
        self.assertEqual(f0 in d, True)
        self.assertEqual(f4 in l, True)
        self.assertEqual(f4 in d, True)
        self.assertEqual(f1 in l, False)
        self.assertEqual(f1 in d, False)

    def testSymbol(self):
        "Test simple symbol instance."

        s0 = Symbol("x", BASIS)
        s1 = Symbol("y", IP)
        s2 = Symbol("z", GEO)
        s3 = Symbol("z", GEO)
        s4 = Symbol("z", IP)

#        print "\nTesting Symbols"
#        print "s0: '%s'" %s0
#        print "s1: '%s'" %s1
#        print "s2: '%s'" %s2
#        print "s3: '%s'" %s3
#        print "s4: '%s'" %s4

        self.assertEqual(repr(s0), "Symbol('x', BASIS)")
        self.assertEqual(repr(s1), "Symbol('y', IP)")
        self.assertEqual(repr(s2), "Symbol('z', GEO)")
        self.assertEqual(repr(s4), "Symbol('z', IP)")

        self.assertEqual(s2 == s3, True)
        self.assertEqual(s2 == s1, False)
        self.assertEqual(s2 == s4, False)
        self.assertEqual(s2 != s3, False)
        self.assertEqual(s2 != s1, True)

        self.assertEqual(s0 < s1, True)
        self.assertEqual(s4 > s1, True)

        self.assertEqual(s0.ops(), 0)
        self.assertEqual(s1.ops(), 0)
        self.assertEqual(s2.ops(), 0)
        self.assertEqual(s3.ops(), 0)
        self.assertEqual(s4.ops(), 0)

        # Test hash
        l = [s0]
        d = {s0:0}
        s5 = Symbol('x', BASIS)

        self.assertEqual(s0 in l, True)
        self.assertEqual(s0 in d, True)
        self.assertEqual(s5 in l, True)
        self.assertEqual(s5 in d, True)

    def testProduct(self):
        "Test simple product instance."

        f0 = FloatValue(-2.0)
        f1 = FloatValue(3.0)
        f2 = FloatValue(0)
        f3 = FloatValue(-1)
        f4 = FloatValue(1)
        f5 = FloatValue(-0.5)
        s0 = Symbol("x", BASIS)
        s1 = Symbol("y", GEO)
        s2 = Symbol("z", GEO)

        p0 = Product([])
        p1 = Product([s0])
        p2 = Product([s0, s1])
        p3 = Product([f1, s0, s1])
        p4 = Product([s0, f2, s2])
        p5 = Product([s0, f0, s1, f1, s2])
        p6 = Product([s0, f3, s1])
        p7 = Product([s0, f4, s1])
        p8 = Product([s0, f0, s2, f5])
        p9 = Product([s0, s1])

#        print "\nTesting Products"
#        print "\np0: [] '%s'" % (p0)
#        print "\np1: %s '%s'" % (s0, p1)
#        print "\np2: %s * %s = '%s'" % (s0, s1, p2)
#        print "\np3: %s * %s * %s = '%s'" % (f1, s0, s1, p3)
#        print "\np4: %s * %s * %s = '%s'" % (s0, f2, s2, p4)
#        print "\np5: %s * %s * %s * %s * %s = '%s'" % (s0, f0, s1, f1, s2, p5)
#        print "\np6: %s * %s * %s = '%s'" % (s0, f3, s1, p6)
#        print "\np7: %s * %s * %s = '%s'" % (s0, f4, s1, p7)
#        print "\np8: %s * %s * %s * %s = '%s'" % (s0, f0, s2, f5, p8)

        self.assertEqual(repr(p0), "Product([FloatValue(0)])")
        self.assertEqual(repr(p1), "Product([Symbol('x', BASIS)])")
        self.assertEqual(repr(p3), "Product([FloatValue(3), Symbol('x', BASIS), Symbol('y', GEO)])")
        self.assertEqual(repr(p6), "Product([FloatValue(-1), Symbol('x', BASIS), Symbol('y', GEO)])")
        self.assertEqual(repr(p7), "Product([Symbol('x', BASIS), Symbol('y', GEO)])")
        self.assertEqual(repr(p8), "Product([Symbol('x', BASIS), Symbol('z', GEO)])")
        self.assertEqual(str(p2), 'x*y')
        self.assertEqual(str(p4), '0')
        self.assertEqual(str(p5), '-6*x*y*z')
        self.assertEqual(str(p6), '-x*y')
        self.assertEqual(str(p7), 'x*y')
        self.assertEqual(str(p8), 'x*z')

        self.assertEqual(p1 == p1, True)
        self.assertEqual(p1 == p7, False)
        self.assertEqual(p4 != p3, True)
        self.assertEqual(p2 == p9, True)
        self.assertEqual(p2 == p3, False)

        self.assertEqual(p0.ops(), 0)
        self.assertEqual(p1.ops(), 0)
        self.assertEqual(p2.ops(), 1)
        self.assertEqual(p3.ops(), 2)
        self.assertEqual(p4.ops(), 0)
        self.assertEqual(p5.ops(), 3)
        self.assertEqual(p6.ops(), 1)
        self.assertEqual(p7.ops(), 1)
        self.assertEqual(p8.ops(), 1)

        # Test hash
        l = [p3]
        d = {p3:0}
        p10 = Product([f1, s0, s1])

        self.assertEqual(p3 in l, True)
        self.assertEqual(p3 in d, True)
        self.assertEqual(p10 in l, True)
        self.assertEqual(p10 in d, True)
        self.assertEqual(p2 in l, False)
        self.assertEqual(p2 in d, False)


    def testSum(self):
        "Test simple sum instance."

        f0 = FloatValue(-2.0)
        f1 = FloatValue(3.0)
        f2 = FloatValue(0)
        s0 = Symbol("x", BASIS)
        s1 = Symbol("y", GEO)
        s2 = Symbol("z", GEO)

        S0 = Sum([])
        S1 = Sum([s0])
        S2 = Sum([s0, s1])
        S3 = Sum([s0, s0])
        S4 = Sum([f0, s0])
        S5 = Sum([s0, f0, s0])
        S6 = Sum([s0, f0, s0, f1])
        S7 = Sum([s0, f0, s1, f2])
        S8 = Sum([s0, f1, s0])
        S9 = Sum([f0, f0, f0, f1, f1, s1])
        S10 = Sum([s1, s0])

#        print "\nTesting Sum"
#        print "\nS0: [] '%s'" % (S0)
#        print "\nS1: %s = '%s'" %(s0, S1)
#        print "\nS2: %s + %s  = '%s'" %(s0, s1, S2)
#        print "\nS3: %s + %s  = '%s'" %(s0, s0, S3)
#        print "\nS4: %s + %s  = '%s'" %(f0, s0, S4)
#        print "\nS5: %s + %s + %s = '%s'" %(s0, f0, s0, S5)
#        print "\nS6: %s + %s + %s + %s = '%s'" %(s0, f0, s0, f1, S6)
#        print "\nS7: %s + %s + %s + %s = '%s'" %(s0, f0, s1, f2, S7)
#        print "\nS8: %s + %s + %s = '%s'" %(s0, f1, s0, S8)
#        print "\nS9: %s + %s + %s + %s + %s + %s = '%s'" %(f0, f0, f0, f1, f1, s1, S9)
#        print "\nS10: %s + %s  = '%s'" %(s1, s0, S10)

        self.assertEqual(repr(S0), "Sum([FloatValue(0)])")
        self.assertEqual(S0.t, CONST)
        self.assertEqual(repr(S1), "Sum([Symbol('x', BASIS)])")
        self.assertEqual(repr(S4), "Sum([Symbol('x', BASIS), FloatValue(-2)])")
        self.assertEqual(repr(S9), "Sum([Symbol('y', GEO)])")

        self.assertEqual(str(S2), "(x + y)")
        self.assertEqual(str(S3), "(x + x)")
        self.assertEqual(str(S5), "(x + x-2)")
        self.assertEqual(str(S6), "(1 + x + x)")
        self.assertEqual(str(S7), "(x + y-2)")
        self.assertEqual(str(S8), "(3 + x + x)")
        self.assertEqual(str(S9), "y")
 
        self.assertEqual(S2 == S2, True)
        self.assertEqual(S2 == S3, False)
        self.assertEqual(S5 != S6, True)
        self.assertEqual(S2 == S10, True)

        self.assertEqual(S0.ops(), 0)
        self.assertEqual(S1.ops(), 0)
        self.assertEqual(S2.ops(), 1)
        self.assertEqual(S3.ops(), 1)
        self.assertEqual(S4.ops(), 1)
        self.assertEqual(S5.ops(), 2)
        self.assertEqual(S6.ops(), 2)
        self.assertEqual(S7.ops(), 2)
        self.assertEqual(S8.ops(), 2)
        self.assertEqual(S9.ops(), 0)

        # Test hash
        l = [S2]
        d = {S2:0}

        self.assertEqual(S2 in l, True)
        self.assertEqual(S2 in d, True)
        self.assertEqual(S10 in l, True)
        self.assertEqual(S10 in d, True)

    def testFraction(self):
        "Test simple fraction instance."

        f0 = FloatValue(-2.0)
        f1 = FloatValue(3.0)
        f2 = FloatValue(0)
        s0 = Symbol("x", BASIS)
        s1 = Symbol("y", GEO)

        F0 = Fraction(f1, f0)
        F1 = Fraction(f2, f0)
        F2 = Fraction(s0, s1)
        F3 = Fraction(s0, f1)
        F4 = Fraction(f0, s1)
        F5 = Fraction(f2, s1)
        F6 = Fraction(s0, s1)

#        print "\nTesting Fractions"
#        print "F0 = frac(%s, %s) = '%s'" %(f1, f0, F0)
#        print "F1 = frac(%s, %s) = '%s'" %(f2, f0, F1)
#        print "F2 = frac(%s, %s) = '%s'" %(s0, s1, F2)
#        print "F3 = frac(%s, %s) = '%s'" %(s0, f1, F3)
#        print "F4 = frac(%s, %s) = '%s'" %(f0, s1, F4)
#        print "F5 = frac(%s, %s) = '%s'" %(f2, s1, F5)
#        print "F6 = frac(%s, %s) = '%s'" %(s0, s1, F6)

        self.assertRaises(RuntimeError, Fraction, f0, f2)
        self.assertRaises(RuntimeError, Fraction, s0, f2)

        self.assertEqual(repr(F0), "Fraction(FloatValue(-1.5), FloatValue(1))")
        self.assertEqual(repr(F2), "Fraction(Symbol('x', BASIS), Symbol('y', GEO))")

        self.assertEqual(str(F0), "-1.5")
        self.assertEqual(str(F1), "0")
        self.assertEqual(str(F2), "x/y")
        self.assertEqual(str(F3), "x/3")
        self.assertEqual(str(F4), "-2/y")
        self.assertEqual(str(F5), "0")

        self.assertEqual(F2 == F2, True)
        self.assertEqual(F2 == F3, False)
        self.assertEqual(F5 != F4, True)
        self.assertEqual(F2 == F6, True)

        self.assertEqual(F0.ops(), 0)
        self.assertEqual(F1.ops(), 0)
        self.assertEqual(F2.ops(), 1)
        self.assertEqual(F3.ops(), 1)
        self.assertEqual(F4.ops(), 1)
        self.assertEqual(F5.ops(), 0)

        # Test hash
        l = [F2]
        d = {F2:0}

        self.assertEqual(F2 in l, True)
        self.assertEqual(F2 in d, True)
        self.assertEqual(F6 in l, True)
        self.assertEqual(F6 in d, True)

    def testMixedSymbols(self):

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
        self.assertEqual(str(mpp1), '3*x*x*y*y')
        self.assertEqual(str(mpp2), '0')
        self.assertEqual(str(mpp3), '9*x*x*x*y*y*y')
        self.assertEqual(str(mps0), 'x*(x + y)')
        self.assertEqual(str(mps1), '(x + x)*(x + y)')
        self.assertEqual(str(mps2), '(x-2)*(x + x-2)')
        self.assertEqual(str(mps3), '(x + x)*(x + x)*(x + y)')
        self.assertEqual(str(mpf0), 'x*x/y')
        self.assertEqual(str(mpf1), 'x/3*x/y')
        self.assertEqual(str(mpf2), '-2/y*x/3')
        self.assertEqual(str(mpf3), 'x/3*x/y*x/y')


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
        self.assertEqual(str(msp4), '0')

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
        self.assertEqual(str(msp1), '(3*x*y + x*y)')
        self.assertEqual(str(msp2), '-6*x*y*z')
        self.assertEqual(str(msp3), '(3*x*y + 3*x*y + x*y)')
        self.assertEqual(str(mss0), '(x + x + y)')
        self.assertEqual(str(mss1), '(x + x + x + y)')
        self.assertEqual(str(mss2), '(x + x + x-4)')
        self.assertEqual(str(mss3), '(x + x + (x + x)*(x + y))')
        self.assertEqual(str(msf0), '(x + x/y)')
        self.assertEqual(str(msf1), '(x/3 + x/y)')
        self.assertEqual(str(msf2), '(x/3-2/y)')
        self.assertEqual(str(msf3), '(x/3 + x/y + x/y)')


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
        self.assertEqual(str(mfp1), '3*x*y/(x*y)')
        self.assertEqual(str(mfp2), '0')
        self.assertEqual(str(mfp3), '3*x*y/(3*x*y/(x*y))')
        self.assertEqual(str(mfs0), '(x + y)/x')
        self.assertEqual(str(mfs1), '(x + x)/(x + y)')
        self.assertEqual(str(mfs2), '(x-2)/(x + x-2)')
        self.assertEqual(str(mfs3), '(x + x)/((x + x)/(x + y))')
        self.assertEqual(str(mff0), '(x/y)/x')
        self.assertEqual(str(mff1), '(x/y)/(x/3)')
        self.assertEqual(str(mff2), '(x/3)/(-2/y)')
        self.assertEqual(str(mff3), '(x/y)/((x/y)/(x/3))')

        # Use p1 as a base expression for Symbol
        s3 = Symbol(get_format()["cos"](p1), CONST, p1, 1)
        self.assertEqual(str(s3), 'std::cos(3*x*y)')
        self.assertEqual(s3.ops(), 3)

    def testOperators(self):
        "Test binary operators"

        f0 = FloatValue(2.0)
        f1 = FloatValue(-3.0)

        s0 = Symbol("x", GEO)
        s1 = Symbol("y", GEO)
        s2 = Symbol("z", GEO)

        p0 = Product([f0, s0])
        p1 = Product([s0, s1])
        p2 = Product([f0, s2])

        S0 = Sum([s0, s1])
        S1 = Sum([s0, s2])

        F0 = Fraction(f0, s1)
        F1 = Fraction(s0, s1)
        F2 = Fraction(s0, S0)
        F3 = Fraction(s0, s1)

        # Test FloatValue '+', only defined for addition by other floats
        # or if self.val == 0
        self.assertEqual(str(f0+f1), '-1')
        self.assertEqual(str(f0+f1+f1+f0+f0), '0')
        self.assertEqual(str(FloatValue(0)+p0), str(p0))
        self.assertRaises(RuntimeError, f1.__add__, p0)

        # Test FloatValue '*', only need one because all other cases are
        # handled by 'other'
        self.assertEqual(str(f0*f0), '4')

        # Test FloatValue '/'
        self.assertEqual(str(f1/f0), '-1.5')
        self.assertEqual(str(f0/s1), '2/y')
        self.assertEqual(str(f1/p0), '-1.5/x')
        self.assertEqual(str(f1/S0), '-3/(x + y)')
        self.assertRaises(RuntimeError, f1.__div__, F0)
        self.assertRaises(RuntimeError, f1.__div__, FloatValue(0))
        self.assertRaises(RuntimeError, f1.__div__, Product([FloatValue(0), s1]))

        # Test Symbol '+', only supported for two equal symbols x + x -> 2*x
        self.assertEqual(str(s0+s0), '2*x')
        self.assertRaises(RuntimeError, s0.__add__, s1)
        self.assertRaises(RuntimeError, s0.__add__, f0)

        # Test Symbol '*', only need to test float, symbol and product. Sum and
        # fraction are handled by 'other'
        self.assertEqual(str(s0*f0), '2*x')
        self.assertEqual(str(s0*s1), 'x*y')
        self.assertEqual(str(s0*p1), 'x*x*y')

        # Test Symbol '/'
        self.assertEqual(str(s0/f0), '0.5*x')
        self.assertEqual(str(s0/s0), '1')
        self.assertEqual(str(s0/s1), 'x/y')
        self.assertEqual(str(s0/S0), 'x/(x + y)')
        self.assertEqual(str(s0/p0), '0.5')
        self.assertEqual(str(s1/p1), '1/x')
        self.assertEqual(str(s2/p0), '0.5*z/x')
        self.assertEqual(str(s2/p1), 'z/(x*y)')
        self.assertRaises(RuntimeError, s0.__div__, F0)
        self.assertRaises(RuntimeError, s1.__div__, FloatValue(0))

        # Test Product '+', only supported for products e.g.,
        # 2*x*y + -5*x*y -> -3*x*y. (and symbols, but only
        # for x + 2*x -> 3*x
        self.assertEqual(str(p0+s0), '3*x')
        self.assertEqual(str(p0+p0), '4*x')
        self.assertEqual(p0+Product([FloatValue(-1), s0]), Symbol('x', GEO))
        self.assertEqual(Product([FloatValue(-1), s0])+s0, FloatValue(0))
        self.assertEqual(str(s0+Product([FloatValue(-1), s0])), '0')
        self.assertRaises(RuntimeError, p0.__add__, FloatValue(0))
        self.assertRaises(RuntimeError, p0.__add__, s2)

        # Test Product '*', only need to test float, symbol and product.
        # Sum and fraction are handled by 'other'
        self.assertEqual(str(p0*FloatValue(0)), '0')
        self.assertEqual(str(p0*f1), '-6*x')
        self.assertEqual(str(p0*s1), '2*x*y')
        self.assertEqual(str(p0*p1), '2*x*x*y')

        # Test Product '/'
        self.assertEqual(str(Product([FloatValue(0), s0])/s0), '0')
        self.assertEqual(str(p0/S0), '2*x/(x + y)')
        self.assertEqual(p1/s1, s0)
        self.assertEqual(p1/p2, Fraction(Product([p1, FloatValue(0.5)]), s2))
        self.assertEqual(p1/s2, Fraction(p1, s2))
        self.assertEqual(p0/Product([f0, p1]), Fraction(FloatValue(1), s1))
        self.assertEqual(p1/p0, Product([FloatValue(0.5), s1]))
        self.assertEqual(p1/p1, FloatValue(1))
        self.assertRaises(RuntimeError, p0.__div__, FloatValue(0))
        self.assertRaises(RuntimeError, p0.__div__, F0)

        # Test Sum '*'
        self.assertEqual(str(S0*f0), '(2*x + 2*y)')
        self.assertEqual(str(S0*s0), '(x*x + x*y)')
        self.assertEqual(str(S0*p0), '(2*x*x + 2*x*y)')
        self.assertEqual(str(S0*S0), '(2*x*y + x*x + y*y)')
        self.assertEqual(str(S0*F0), '(2 + 2*x/y)')

        # Test Sum '/'
        self.assertEqual(str(S0/f0), '(0.5*x + 0.5*y)')
        self.assertEqual(str(S0/s0), '(1 + y/x)')
        self.assertEqual(str(S0/p0), '(0.5 + 0.5*y/x)')
        self.assertEqual(str(S0/p1), '(1/x + 1/y)')
        self.assertEqual(str(S0/S0), '(x + y)/(x + y)')
        self.assertEqual(str(S0/S1), '(x + y)/(x + z)')
        self.assertRaises(RuntimeError, S0.__div__, FloatValue(0))
        self.assertRaises(RuntimeError, S0.__div__, F0)

        # Test Fraction '+'
        self.assertEqual(str(F1+F3), '2*x/y')
        self.assertEqual(str(F0+F1), '(2 + x)/y')
        self.assertRaises(RuntimeError, F2.__add__, FloatValue(0))
        self.assertRaises(RuntimeError, F2.__add__, F3)

        # Test Fraction '*'
        self.assertEqual(str(F1*f0), '2*x/y')
        self.assertEqual(str(F1*s0), 'x*x/y')
        self.assertEqual(str(F1*p1), 'x*x')
        self.assertEqual(str(F1*S0), '(x + x*x/y)')
        self.assertEqual(repr(F1*S0), repr(Sum([s0, Fraction(Product([Symbol('x', GEO), Symbol('x', GEO)]), Symbol('y', GEO))])))
        self.assertEqual(str(F1*F0), '2*x/(y*y)')


    def testExpandOperations(self):
        s0 = Product([FloatValue(-1), Symbol("x", GEO)])
        s1 = Symbol("y", GEO)
        s2 = Product([FloatValue(5), Symbol("z", IP)])
        s3 = Product([FloatValue(-4), Symbol("z", GEO)])
        f0 = FloatValue(-1)
        f1 = FloatValue(2)
        f2 = FloatValue(1)

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

        F4x = F4.expand()
        F5x = F5.expand()
        F6x = F6.expand()
        F7x = F7.expand()

#        print "\nF4: '%s'" %F4
#        print "F4x: '%s'" %F4x
#        print "\nF5: '%s'" %F5
#        print "F5x: '%s'" %F5x
#        print "\nF6: '%s'" %F6
#        print "F6x: '%s'" %F6x
#        print "\nF7: '%s'" %F7
#        print "F7x: '%s'" %F7x

        self.assertAlmostEqual(eval(str(F4)), eval(str(F4x)))
        self.assertAlmostEqual(eval(str(F5)), eval(str(F5x)))
        self.assertAlmostEqual(eval(str(F6)), eval(str(F6x)))
        self.assertAlmostEqual(eval(str(F7)), eval(str(F7x)))

        self.assertEqual(F4.ops(), 5)
        self.assertEqual(F4x.ops(), 1)
        self.assertEqual(F5.ops(), 6)
        self.assertEqual(F5x.ops(), 5)
        self.assertEqual(F6.ops(), 9)
        self.assertEqual(F6x.ops(), 1)
        self.assertEqual(F7.ops(), 2)
        self.assertEqual(F7x.ops(), 1)

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

        P4 = Product([s1, Sum([s0, s1])])
        P5 = Product([s0, E0])
        P6 = Product([s1])
        S4 = Sum([s1])


        # Create 'real' term that caused me trouble
        P00 = Product([Symbol("Jinv_00", GEO)]*2)
        P01 = Product([Symbol("Jinv_01", GEO)]*2)
        P20 = Product([Symbol("Jinv_00", GEO), Product([f1, Symbol("Jinv_20", GEO)]) ])
        P21 = Product([Symbol("Jinv_01", GEO), Product([f1, Symbol("Jinv_21", GEO)]) ])
        PS0 = Product([Symbol("Jinv_22", GEO), Sum([P00, P01])])
        PS1 = Product([ Product([f0, Symbol("Jinv_02", GEO)]), Sum([P20, P21])])
        SP = Sum([PS0, PS1])

        PFx = PF.expand()
        E0x = E0.expand()
        E1x = E1.expand()
        E2x = E2.expand()
        E3x = E3.expand()
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
        self.assertEqual(SP.ops(), 11)
        self.assertEqual(SPx.ops(), 13)
        self.assertEqual(P4.ops(), 2)
        self.assertEqual(P4x.ops(), 3)
        self.assertEqual(P5.ops(), 5)
        self.assertEqual(P5x.ops(), 9)

    def testReduceVarType(self):

        f1 = FloatValue(1)
        f2 = FloatValue(2)
        f3 = FloatValue(3)
        f5 = FloatValue(5)
        fm4 = FloatValue(-4)

        B0 = Symbol("B0",BASIS)
        B1 = Symbol("B1", BASIS)
        Bm4 = Product([fm4, B1])
        B5 = Product([f5, B0])

        I0 = Symbol("I0", IP)
        I5 = Product([f5, I0])

        G0 = Symbol("G0", GEO)
        G3 = Product([f3, G0])


        C0 = Symbol("C0", CONST)
        C2 = Product([f2, C0])

        p0 = Product([B0,I5])
        p1 = Product([B0,B1])

        S0 = Sum([B0, I5])
        S1 = Sum([p0, p1])
        S2 = Sum([B0, B1])
        S3 = Sum([B0, p0])
        S4 = Sum([f5, p0])

        F0 = Fraction(B0,I5).expand()
        F1 = Fraction(p1,I5).expand()
        F2 = Fraction(G3,S2).expand()
        F3 = Fraction(G3,S3).expand()

        r0 = B0.reduce_vartype(BASIS)
        r1 = B0.reduce_vartype(CONST)

        rp0 = p0.reduce_vartype(BASIS)
        rp1 = p0.reduce_vartype(IP)
        rp2 = p1.reduce_vartype(BASIS)
        rp3 = p1.reduce_vartype(GEO)

        rs0 = S0.reduce_vartype(BASIS)
        rs1 = S0.reduce_vartype(IP)
        rs2 = S1.reduce_vartype(BASIS)
        rs3 = S4.reduce_vartype(BASIS)
        rs4 = S4.reduce_vartype(CONST)

        rf0 = F0.reduce_vartype(BASIS)
        rf1 = F1.reduce_vartype(BASIS)
        rf2 = F0.reduce_vartype(IP)
        rf3 = F2.reduce_vartype(BASIS)
        rf4 = F3.reduce_vartype(BASIS)


#        print "%s, red(BASIS): ('%s', '%s')" %(B0, r0[0], r0[1])
#        print "%s, red(CONST): ('%s', '%s')" %(B0, r1[0], r1[1])

#        print "\n%s, red(BASIS): ('%s', '%s')" %(p0, rp0[0], rp0[1])
#        print "%s, red(IP):    ('%s', '%s')" %(p0, rp1[0], rp1[1])
#        print "%s, red(BASIS): ('%s', '%s')" %(p1, rp2[0], rp2[1])
#        print "%s, red(CONST): ('%s', '%s')" %(p1, rp3[0], rp3[1])

#        print "\n%s, red(BASIS): ('%s', '%s')" %(S0, rs0[0], rs0[1])
#        print "%s, red(IP):    ('%s', '%s')" %(S0, rs1[0], rs1[1])
#        print "%s, red(BASIS): '%s', '%s'" %(S1, rs2[0], rs2[1])
#        print "%s, red(BASIS): '%s', '%s'" %(S4, rs3[0], rs3[1])
#        print "%s, red(BASIS): '%s'" %(S4, rs4[0])

#        print "\n%s, red(BASIS): ('%s', '%s')" %(F0, rf0[0], rf0[1])
#        print "%s, red(BASIS): ('%s', '%s')" %(F1, rf1[0], rf1[1])
#        print "%s, red(IP): ('%s', '%s')" %(F0, rf2[0], rf2[1])
#        print "%s, red(BASIS): ('%s', '%s')" %(F2, rf3[0], rf3[1])
#        print "%s, red(BASIS): ('%s', '%s')" %(F3, rf4[0], rf4[1])

        self.assertEqual((B0, f1), r0)
        self.assertEqual(((), B0), r1)

        self.assertEqual((B0, I5), rp0)
        self.assertEqual((I0, B5),  rp1)
        self.assertEqual((p1, f1), rp2)
        self.assertEqual(((), p1),  rp3)

        self.assertEqual(((), I5), rs0[0])
        self.assertEqual((B0, f1), rs0[1])
        self.assertEqual((I0, f5), rs1[0])
        self.assertEqual(((), B0), rs1[1])
        self.assertEqual((Product([B0, B1]), f1), rs2[0])
        self.assertEqual((B0, I5), rs2[1])
        self.assertEqual(((), f5), rs3[0])
        self.assertEqual((B0, I5), rs3[1])
        self.assertEqual((f5, Sum([f1, Product([B0, I0])])), rs4[0])

        self.assertEqual((B0, Fraction(FloatValue(0.2), I0)), rf0)
        self.assertEqual((Product([B0, B1]), Fraction(FloatValue(0.2), I0)), rf1)
        self.assertEqual( ( Fraction(f1, I0), Product([FloatValue(0.2), B0]) ), rf2)
        self.assertEqual(((), F2), rf3)
        self.assertEqual( ( Fraction(f1, B0), Fraction( G3, Sum([I5, f1]))), rf4)

    def testReduceOperations(self):

        print "\nTesting ReduceOperations"

        # Aux. variables
        f2 = FloatValue(2)
        f0_5 = FloatValue(0.5)
        f1 = FloatValue(1.0)
        fm1 = FloatValue(-1.0)

        x = Symbol("x", GEO)
        y = Symbol("y", GEO)
        z = Symbol("z", GEO)
        a = Symbol("a", GEO)
        b = Symbol("b", GEO)
        c = Symbol("c", GEO)
        d = Symbol("d", GEO)

        # Simple expand and reduce simple float and symbol objects
        fx2 = f2.expand()
        xx = x.expand()

        fr2 = fx2.reduce_ops()
        xr = xx.reduce_ops()

#        print "\nTest float and symbol"
#        print "f0:  '%s'" %f2
#        print "fx0: '%s'" %fx2
#        print "fr0: '%s'" %fr2
#        print
#        print "x:  '%s'" %x
#        print "xx: '%s'" %xx
#        print "xr: '%s'" %xr

        self.assertEqual(f2, fr2)
        self.assertEqual(x, xr)

        # Test product
        p0 = f2*x
        p1 = y*x
        p2 = x*f2/y
        p3 = x*Sum([x, y])

        px0 = p0.expand()
        px1 = p1.expand()

        pr0 = px0.reduce_ops()
        pr1 = px1.reduce_ops()

#        print "\nTest product"
#        print "p0:  '%s'" %p0
#        print "px0: '%s'" %px0
#        print "pr0: '%s'" %pr0
#        print
#        print "p1:  '%s'" %p1
#        print "px1: '%s'" %px1
#        print "pr1: '%s'" %pr1

        self.assertEqual(p0, pr0)
        self.assertEqual(p1, pr1)

        # Test fraction
        F0 = Fraction(p0, y)
        F1 = Fraction(x, p0)
        F2 = Fraction(p0, p1)
        F3 = Fraction(Sum([x*x, x*y]), y)
        F4 = Fraction(Sum([f2*x, x*y]), a)

        Fx0 = F0.expand()
        Fx1 = F1.expand()
        Fx2 = F2.expand()
        Fx3 = F3.expand()
        Fx4 = F4.expand()

        Fr0 = Fx0.reduce_ops()
        Fr1 = Fx1.reduce_ops()
        Fr2 = Fx2.reduce_ops()
        Fr3 = Fx3.reduce_ops()
        Fr4 = Fx4.reduce_ops()

#        print "\nTest fraction"
#        print "F0:  '%s'" %F0
#        print "Fx0: '%s'" %Fx0
#        print "Fr0: '%s'" %Fr0
#        print
#        print "F1:  '%s'" %F1
#        print "Fx1: '%s'" %Fx1
#        print "Fr1: '%s'" %Fr1
#        print
#        print "F2:  '%s'" %F2
#        print "Fx2: '%s'" %Fx2
#        print "Fr2: '%s'" %Fr2
#        print
#        print "F3:  '%s'" %F3
#        print "Fx3: '%s'" %Fx3
#        print "Fr3: '%s'" %Fr3
#        print
#        print "F4:  '%s'" %F4
#        print "Fx4: '%s'" %Fx4
#        print "Fr4: '%s'" %Fr4

        self.assertEqual(Fr0, F0)
        self.assertEqual(Fr1, f0_5)
        self.assertEqual(Fr2, Fraction(f2, y))
        self.assertEqual(str(Fr3), "x*(1 + x/y)")
        self.assertEqual(str(Fr4), "x*(2 + y)/a")

        # Test sum
        # TODO: Here we might have to add additional tests
        S0 = Sum([x, y])
        S1 = Sum([p0, p1])
        S2 = Sum([x, p1])
        S3 = Sum([p0, f2*y])
        S4 = Sum([f2*p1, z*p1])
        S5 = Sum([x, x*x, x*x*x])
        S6 = Sum([a*x*x, b*x*x*x, c*x*x, d*x*x*x])
        S7 = Sum([p0, p1, x*x, f2*z, y*z])
        S8 = Sum([a*y, b*y, x*x*x*y, x*x*x*z])
        S9 = Sum([a*y, b*y, c*y, x*x*x*y, f2*x*x, x*x*x*z])
        S10 = Sum([f2*x*x*y, x*x*y*z])
        S11 = Sum([f2*x*x*y*y, x*x*y*y*z])
        S12 = Sum([f2*x*x*y*y, x*x*y*y*z, a*z, b*z, c*z])
        S13 = Sum([Fraction(f1, x), Fraction(f1, y)])
        S14 = Sum([Fraction(fm1, x), Fraction(fm1, y)])
        S15 = Sum([Fraction(f2, x), Fraction(f2, x)])
        S16 = Sum([Fraction(f2*x, y*z), Fraction(f0_5, y*z)])
        S17 = Sum([(f2*x*y)/a, (x*y*z)/b])
        S18 = Sum([(x*y)/a, (x*z)/a, f2/a, (f2*x*y)/a])
        S19 = Sum([(f2*x)/a, (x*y)/a, z*x])
        S20 = Product([ Sum([x, y]), Fraction(a, b), Fraction( Product([c, d]), z ) ])

        Sx0 = S0.expand()
        Sx1 = S1.expand()
        Sx2 = S2.expand()
        Sx3 = S3.expand()
        Sx4 = S4.expand()
        Sx5 = S5.expand()
        Sx6 = S6.expand()
        Sx7 = S7.expand()
        Sx8 = S8.expand()
        Sx9 = S9.expand()
        Sx10 = S10.expand()
        Sx11 = S11.expand()
        Sx12 = S12.expand()
        Sx13 = S13.expand()
        Sx14 = S14.expand()
        Sx15 = S15.expand()
        Sx16 = S16.expand()
        Sx17 = S17.expand()
        Sx18 = S18.expand()
        Sx19 = S19.expand()
        Sx20 = S20.expand()

        Sr0 = Sx0.reduce_ops()
        Sr1 = Sx1.reduce_ops()
        Sr2 = Sx2.reduce_ops()
        Sr3 = Sx3.reduce_ops()
        Sr4 = Sx4.reduce_ops()
        Sr5 = Sx5.reduce_ops()
        Sr6 = Sx6.reduce_ops()
        Sr7 = Sx7.reduce_ops()
        Sr8 = Sx8.reduce_ops()
        Sr9 = Sx9.reduce_ops()
        Sr10 = Sx10.reduce_ops()
        Sr11 = Sx11.reduce_ops()
        Sr12 = Sx12.reduce_ops()
        Sr13 = Sx13.reduce_ops()
        Sr14 = Sx14.reduce_ops()
        Sr15 = Sx15.reduce_ops()
        Sr16 = Sx16.reduce_ops()
        Sr17 = Sx17.reduce_ops()
        Sr18 = Sx18.reduce_ops()
        Sr19 = Sx19.reduce_ops()
        Sr20 = Sx20.reduce_ops()

#        print "Test sum"
#        print "S0:  '%s'" %S0
#        print "Sx0: '%s'" %Sx0
#        print "Sr0: '%s'" %Sr0
#        print
#        print "S1:  '%s'" %S1
#        print "Sx1: '%s'" %Sx1
#        print "Sr1: '%s'" %Sr1
#        print
#        print "S2:  '%s'" %S2
#        print "Sx2: '%s'" %Sx2
#        print "Sr2: '%s'" %Sr2
#        print
#        print "S3:  '%s'" %S3
#        print "Sx3: '%s'" %Sx3
#        print "Sr3: '%s'" %Sr3
#        print
#        print "S4:  '%s'" %S4
#        print "Sx4: '%s'" %Sx4
#        print "Sr4: '%s'" %Sr4
#        print
#        print "S5:  '%s'" %S5
#        print "Sx5: '%s'" %Sx5
#        print "Sr5: '%s'" %Sr5
#        print
#        print "S6:  '%s'" %S6
#        print "Sx6: '%s'" %Sx6
#        print "Sr6: '%s'" %Sr6
#        print
#        print "S7:  '%s'" %S7
#        print "Sx7: '%s'" %Sx7
#        print "Sr7: '%s'" %Sr7
#        print
#        print "S8:  '%s'" %S8
#        print "Sx8: '%s'" %Sx8
#        print "Sr8: '%s'" %Sr8
#        print
#        print "S9:  '%s'" %S9
#        print "Sx9: '%s'" %Sx9
#        print "Sr9: '%s'" %Sr9
#        print
#        print "S10:  '%s'" %S10
#        print "Sx10: '%s'" %Sx10
#        print "Sr10: '%s'" %Sr10
#        print
#        print "S11:  '%s'" %S11
#        print "Sx11: '%s'" %Sx11
#        print "Sr11: '%s'" %Sr11
#        print
#        print "S12:  '%s'" %S12
#        print "Sx12: '%s'" %Sx12
#        print "Sr12: '%s'" %Sr12
#        print
#        print "S13:  '%s'" %S13
#        print "Sx13: '%s'" %Sx13
#        print "Sr13: '%s'" %Sr13
#        print
#        print "S14:  '%s'" %S14
#        print "Sx14: '%s'" %Sx14
#        print "Sr14: '%s'" %Sr14
#        print
#        print "S15:  '%s'" %S15
#        print "Sx15: '%s'" %Sx15
#        print "Sr15: '%s'" %Sr15
#        print
#        print "S16:  '%s'" %S16
#        print "Sx16: '%s'" %Sx16
#        print "Sr16: '%s'" %Sr16
#        print
#        print "S17:  '%s'" %S17
#        print "Sx17: '%s'" %Sx17
#        print "Sr17: '%s'" %Sr17
#        print
#        print "S18:  '%s'" %S18
#        print "Sx18: '%s'" %Sx18
#        print "Sr18: '%s'" %Sr18
#        print
#        print "S19:  '%s'" %S19
#        print "Sx19: '%s'" %Sx19
#        print "Sr19: '%s'" %Sr19
#        print
#        print "S20:  '%s'" %S20
#        print "Sx20: '%s'" %Sx20
#        print "Sr20: '%s'" %Sr20

        self.assertEqual(Sr0, S0)
        self.assertEqual(str(Sr1), "x*(2 + y)")
        # TODO: Should this be (x + x*y)?
        self.assertEqual(str(Sr2), "x*(1 + y)")
#        self.assertEqual(str(Sr2), "(x + x*y)")
        self.assertEqual(str(Sr3), "2*(x + y)")
        self.assertEqual(str(Sr4), "x*y*(2 + z)")
        self.assertEqual(str(Sr5), "x*(1 + x*(1 + x))")
        self.assertEqual(str(Sr6), "x*x*(a + c + x*(b + d))")
        self.assertEqual(str(Sr7), "(x*(2 + x + y) + z*(2 + y))")
        self.assertEqual(str(Sr8), "(x*x*x*(y + z) + y*(a + b))")
        self.assertEqual(str(Sr9), "(x*x*(2 + x*(y + z)) + y*(a + b + c))")
        self.assertEqual(str(Sr10), "x*x*y*(2 + z)")
        self.assertEqual(str(Sr11), "x*x*y*y*(2 + z)")
        self.assertEqual(str(Sr12), "(x*x*y*y*(2 + z) + z*(a + b + c))")
        self.assertEqual(str(Sr13), "(1/x + 1/y)")
        self.assertEqual(str(Sr14), "(-1/x-1/y)")
        self.assertEqual(str(Sr15), "4/x")
        self.assertEqual(str(Sr16), "(0.5 + 2*x)/(y*z)")
        self.assertEqual(str(Sr17), "x*y*(2/a + z/b)")
        self.assertEqual(str(Sr18), "(2 + x*(z + 3*y))/a")
        self.assertEqual(str(Sr19), "x*(z + (2 + y)/a)")
        self.assertEqual(str(Sr20), "a*c*d*(x + y)/(b*z)")

        self.assertEqual(S0.ops(), 1)
        self.assertEqual(Sr0.ops(), 1)
        self.assertEqual(S1.ops(), 3)
        self.assertEqual(Sr1.ops(), 2)
        self.assertEqual(S2.ops(), 2)
        self.assertEqual(Sr2.ops(), 2)
        self.assertEqual(S3.ops(), 3)
        self.assertEqual(Sr3.ops(), 2)
        self.assertEqual(S4.ops(), 5)
        self.assertEqual(Sr4.ops(), 3)
        self.assertEqual(S5.ops(), 5)
        self.assertEqual(Sr5.ops(), 4)
        self.assertEqual(S6.ops(), 13)
        self.assertEqual(Sr6.ops(), 6)
        self.assertEqual(S7.ops(), 9)
        self.assertEqual(Sr7.ops(), 6)
        self.assertEqual(S8.ops(), 11)
        self.assertEqual(Sr8.ops(), 7)
        self.assertEqual(S9.ops(), 16)
        self.assertEqual(Sr9.ops(), 9)
        self.assertEqual(S10.ops(), 7)
        self.assertEqual(Sr10.ops(), 4)
        self.assertEqual(S11.ops(), 9)
        self.assertEqual(Sr11.ops(), 5)
        self.assertEqual(S12.ops(), 15)
        self.assertEqual(Sr12.ops(), 9)
        self.assertEqual(S13.ops(), 3)
        self.assertEqual(Sr13.ops(), 3)
        self.assertEqual(S14.ops(), 3)
        self.assertEqual(Sr14.ops(), 3)
        self.assertEqual(S15.ops(), 3)
        self.assertEqual(Sr15.ops(), 1)
        self.assertEqual(S16.ops(), 6)
        self.assertEqual(Sr16.ops(), 4)
        self.assertEqual(S17.ops(), 7)
        self.assertEqual(Sr17.ops(), 5)
        self.assertEqual(S18.ops(), 11)
        self.assertEqual(Sr18.ops(), 5)
        self.assertEqual(S19.ops(), 7)
        self.assertEqual(Sr19.ops(), 4)

    def testNotFinished(self):
        "Stuff that would be nice to implement."

        f0 = FloatValue(4)
        f1 = FloatValue(2)
        f2 = FloatValue(8)
        s0 = Symbol("x", GEO)
        s1 = Symbol("y", GEO)
        s2 = Symbol("z", GEO)
        a = Symbol("a", GEO)
        b = Symbol("b", GEO)
        c = Symbol("c", GEO)

        # Aux. expressions
        p0 = Product([f1, s0])
        p1 = Product([f2, s1])
        p2 = Product([s0, s1])

        F0 = Fraction(f0, s0)

        S0 = Sum([p0, p1])
        S1 = Sum([s0, p2])
        S2 = Sum([FloatValue(1), s1])
        S3 = Sum([F0, F0])

        # Thing to be implemented
        e0 = f0 / S0
        e1 = s0 / S1
        e2 = S2 / S1
        e3 = group_fractions(S3)
        e4 = Sum([Fraction(f1*s0, a*b*c), Fraction(s0, a*b)]).expand().reduce_ops()

        # Tests that pass the current implementation
        self.assertEqual(str(e0), '4/(2*x + 8*y)')
        self.assertEqual(str(e1), 'x/(x + x*y)')
        self.assertEqual(str(e2), '(1 + y)/(x + x*y)')
        self.assertEqual(str(e3), '8/x')
        self.assertEqual(str(e4), 'x*(1/(a*b) + 2/(a*b*c))')

        # Tests that should pass in future implementations (change NotEqual to Equal)
        self.assertNotEqual(str(e0), '2/(x + 4*y)')
        self.assertNotEqual(str(e1), '1/(1 + y)')
        self.assertNotEqual(str(e2), '1/x')
        self.assertNotEqual(str(e4), 'x*(2/c + 1)/(a*b)')


    def testCache(self):
        raise RuntimeError("Tests not implemented")


    def testDGElastoDyn(self):
        expr = Product([
                       Sum([
                            Symbol("F0", IP),
                            Symbol("F1", IP)
                          ]),
                       Fraction(
                                 Symbol("w4", GEO),
                                 Symbol("w3", GEO)
                                ),
                       Fraction(
                                 Product([
                                          Symbol("w2", GEO),
                                          Symbol("w5", GEO)
                                         ]),
                                 Symbol("w6", GEO)
                                )
                      ])

        print
        start = time.time()
        expr_exp = expr.expand()
        print "DGElastoDyn: time, expand():     ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "DGElastoDyn: time, reduce_ops(): ", time.time() - start

        print "expr.ops():     ", expr.ops()
        print "expr_exp.ops(): ", expr_exp.ops()
        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        F0, F1, w2, w3, w4, w5, w6 = (3.12, -8.1, -45.3, 17.5, 2.2, 5.3, 9.145)
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))

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

        print
        start = time.time()
        expr_exp = expr.expand()
        print "ReduceGIP: time, expand()      ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "ReduceGIP: time, reduce_ops(): ", time.time() - start

        print "expr.ops():     ", expr.ops()
        print "expr_exp.ops(): ", expr_exp.ops()
        print "expr_red.ops(): ", expr_red.ops()

#        print "expr: ", expr
#        print "exp:  ", expr_exp
#        print "red:  ", expr_red

        W9 = 9
        F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20 = [0.123 * i for i in range(1,21)]
        G0, G1, G2, G3, G4, G5, G6, G7, G8, G9 = [2.64 + 1.0/i for i in range(20, 30)]

        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(expr.ops() > expr_red.ops(), True)


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

        print
        start = time.time()
        expr_exp = expr.expand()
        print "Poisson: time, expand():     ", time.time() - start

        start = time.time()
        poisson_exp = expand_operations(poisson, get_format())
        print "Poisson: time, old expand(): ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "Poisson: time, reduce_ops(): ", time.time() - start

        start = time.time()
        poisson_red = reduce_operations(poisson, get_format())
        print "Poisson: time, old reduce(): ", time.time() - start

        ops = operation_count(poisson_exp, get_format())
        print "expr.ops():           ", expr.ops()
        print "Poisson old exp: ops: ", ops
        print "expr_exp.ops():       ", expr_exp.ops()

        ops = operation_count(poisson_red, get_format())
        print "Poisson old red: ops: ", ops
        print "expr_red.ops():       ", expr_red.ops()

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

        print
        start = time.time()
        expr_exp = expr.expand()
        print "Elasticity2D: time, expand():     ", time.time() - start

        start = time.time()
        elasticity_exp = expand_operations(elasticity, get_format())
        print "Elasticity2D: time, old expand(): ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "Elasticity2D: time, reduce_ops(): ", time.time() - start

        start = time.time()
        elasticity_red = reduce_operations(elasticity, get_format())
        print "Elasticity2D: time, old reduce(): ", time.time() - start

        ops = operation_count(elasticity_exp, get_format())
        print "expr.ops():                ", expr.ops()
        print "Elasticity2D old exp: ops: ", ops
        print "expr_exp.ops():            ", expr_exp.ops()

        ops = operation_count(elasticity_red, get_format())
        print "Elasticity2D old red: ops: ", ops
        print "expr_red.ops():            ", expr_red.ops()

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

        print
        start = time.time()
        expr_exp = expr.expand()
        print "ElasticityTerm: time, expand():     ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "ElasticityTerm: time, reduce_ops(): ", time.time() - start

        print "expr.ops():     ", expr.ops()
        print "expr_exp.ops(): ", expr_exp.ops()
        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        # Generate code
        ip_consts = {}
        geo_consts = {}
        trans_set = set()

        start = time.time()
        opt_code = optimise_code(expr, ip_consts, geo_consts, trans_set)
        print "ElasticityTerm, optimise_code(): ", time.time() - start

        det, W1, Jinv_00, Jinv_21, FE0_C2_D001_0_j, FE0_C2_D001_0_k = [0.123 + i for i in range(6)]
        G0 = eval(str(geo_consts.items()[0][0]))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))


    def testElasWeighted(self):
        expr = Product([
                          Symbol('W4', IP),
                          Symbol('det', GEO),
                          Sum([
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', BASIS),
                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                        Symbol('Jinv_00', GEO),
                                        Symbol('Jinv_11', GEO),
                                        Symbol('w1', GEO)
                                        ]),
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', BASIS),
                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                        Symbol('Jinv_01', GEO),
                                        Symbol('Jinv_10', GEO),
                                        Symbol('w0', GEO)
                                        ]),
                              Product([
                                        Symbol('w2', GEO),
                                        Sum([
                                              Product([
                                                      Symbol('FE0_C1_D01_ip_j', BASIS),
                                                      Symbol('FE0_C1_D01_ip_k', BASIS),
                                                      Symbol('Jinv_00', GEO),
                                                      Symbol('Jinv_11', GEO),
                                                      Symbol('w1', GEO)
                                                      ]),
                                              Product([
                                                      Symbol('FE0_C1_D01_ip_j', BASIS),
                                                      Symbol('FE0_C1_D01_ip_k', BASIS),
                                                      Symbol('Jinv_01', GEO),
                                                      Symbol('Jinv_10', GEO),
                                                      Symbol('w0', GEO)
                                                      ])
                                            ])
                                      ])
                              ])
                          ])
                                                       
        print
        start = time.time()
        expr_exp = expr.expand()
        print "ElasWeighted: time, expand():     ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "ElasWeighted: time, reduce_ops(): ", time.time() - start

        print "expr.ops():     ", expr.ops()
        print "expr_exp.ops(): ", expr_exp.ops()
        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        # Generate code
        ip_consts = {}
        geo_consts = {}
        trans_set = set()

        start = time.time()
        opt_code = optimise_code(expr, ip_consts, geo_consts, trans_set)
        print "ElasWeighted, optimise_code(): ", time.time() - start

        det, W4, w0, w1, w2, Jinv_00, Jinv_01, Jinv_11, Jinv_10, FE0_C1_D01_ip_j, FE0_C1_D01_ip_k = [0.123 + i for i in range(11)]
        G0 = eval(str(geo_consts.items()[0][0]))
        Gip0 = eval(str(ip_consts.items()[0][0]))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))


    def testElasWeighted2(self):

        expr = Product([
                        Symbol('W4', IP),
                        Symbol('det', GEO),
                        Sum([
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', BASIS),
                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                        Symbol('Jinv_00', GEO),
                                        Symbol('Jinv_10', GEO),
                                        Symbol('w1', GEO)
                                        ]),
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', BASIS),
                                        Symbol('Jinv_01', GEO),
                                        Sum([
                                              Product([
                                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                                        Symbol('Jinv_11', GEO),
                                                        Symbol('w0', GEO)
                                                        ]),
                                              Product([
                                                        FloatValue(2),
                                                        Symbol('FE0_C1_D01_ip_k', BASIS),
                                                        Symbol('Jinv_11', GEO),
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
                                                    Symbol('Jinv_10', GEO),
                                                    Symbol('w1', GEO)
                                                    ]),
                                            Product([
                                                    Symbol('FE0_C1_D01_ip_j', BASIS),
                                                    Symbol('Jinv_01', GEO),
                                                    Sum([
                                                          Product([
                                                                  Symbol('FE0_C1_D01_ip_k', BASIS),
                                                                  Symbol('Jinv_11', GEO),
                                                                  Symbol('w0', GEO)
                                                                  ]),
                                                          Product([
                                                                  FloatValue(2),
                                                                  Symbol('FE0_C1_D01_ip_k', BASIS),
                                                                  Symbol('Jinv_11', GEO),
                                                                  Symbol('w1', GEO)
                                                                  ])
                                                          ])
                                                    ])
                                            ])
                                        ])
                              ])
                        ])

        print
        start = time.time()
        expr_exp = expr.expand()
        print "ElasWeighted2: time, expand():     ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "ElasWeighted2: time, reduce_ops(): ", time.time() - start

        print "expr.ops():     ", expr.ops()
        print "expr_exp.ops(): ", expr_exp.ops()
        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        # Generate code
        ip_consts = {}
        geo_consts = {}
        trans_set = set()

        start = time.time()
        opt_code = optimise_code(expr, ip_consts, geo_consts, trans_set)
        print "ElasWeighted2, optimise_code(): ", time.time() - start

        det, W4, w0, w1, w2, Jinv_00, Jinv_01, Jinv_11, Jinv_10, FE0_C1_D01_ip_j, FE0_C1_D01_ip_k = [0.123 + i for i in range(11)]
        G0 = eval(str(geo_consts.items()[0][0]))
        Gip0 = eval(str(ip_consts.items()[0][0]))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))


def suite():
    suite = unittest.TestSuite()
    # Classes and member functions
    suite.addTest(Tests('testFloat'))
    suite.addTest(Tests('testSymbol'))
    suite.addTest(Tests('testProduct'))
    suite.addTest(Tests('testSum'))
    suite.addTest(Tests('testFraction'))
    suite.addTest(Tests('testMixedSymbols'))
    suite.addTest(Tests('testOperators'))
    suite.addTest(Tests('testExpandOperations'))
    suite.addTest(Tests('testReduceVarType'))
    suite.addTest(Tests('testReduceOperations'))

    # Misc.
    suite.addTest(Tests('testNotFinished'))
#    suite.addTest(Tests('testCache'))

    # 'Real' expressions (expand and reduce)
    suite.addTest(Tests('testDGElastoDyn'))
    suite.addTest(Tests('testReduceGIP'))
    suite.addTest(Tests('testPoisson'))
    suite.addTest(Tests('testElasticity2D'))

    # 'Real' expressions (generate code)
    suite.addTest(Tests('testElasticityTerm'))
    suite.addTest(Tests('testElasWeighted'))
    suite.addTest(Tests('testElasWeighted2'))

    return suite


if __name__ == "__main__":


    if format == None:
        set_format(Format(FFC_OPTIONS).format)

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(suite())

