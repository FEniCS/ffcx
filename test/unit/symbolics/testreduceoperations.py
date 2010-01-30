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
set_float_formatting(FFC_OPTIONS)

class TestReduceOperations(unittest.TestCase):

    def testReduceOperations(self):

        f_1 = format["float"](1)
        f_2 = format["float"](2)

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
        self.assertEqual(str(Fr3), "x*(%s + x/y)" % f_1)
        self.assertEqual(str(Fr4), "x*(%s + y)/a" % f_2)

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
        S21 = Sum([a*x, b*x, c*x, x*y, x*z, f2*y, a*y, b*y, f2*z, a*z, b*z])
        S22 = Sum([ FloatValue(0.5)*x/y, FloatValue(-0.5)*x/y ])
        S23 = Sum([x*y*z, x*y*y*y*z*z*z, y*y*y*z*z*z*z, z*z*z*z*z])


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
        Sx21 = S21.expand()
        Sx22 = S22.expand()
        Sx23 = S23.expand()

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
        Sr21 = Sx21.reduce_ops()
        Sr22 = Sx22.reduce_ops()
        Sr23 = Sx23.reduce_ops()

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
#        print
#        print "S21:  '%s'" %S21
#        print "Sx21: '%s'" %Sx21
#        print "Sr21: '%s'" %Sr21
#        print
#        print "S22:  '%s'" %S22
#        print "Sx22: '%s'" %Sx22
#        print "Sr22: '%s'" %Sr22
#        print
#        print "S23:  '%s'" %S23
#        print "Sx23: '%s'" %Sx23
#        print "Sr23: '%s'" %Sr23

        self.assertEqual(Sr0, S0)
        self.assertEqual(str(Sr1), "x*(%s + y)" % f_2)
        # TODO: Should this be (x + x*y)?
        self.assertEqual(str(Sr2), "x*(%s + y)" % f_1)
#        self.assertEqual(str(Sr2), "(x + x*y)")
        self.assertEqual(str(Sr3), "%s*(x + y)" % f_2)
        self.assertEqual(str(Sr4), "x*y*(%s + z)" % f_2)
        self.assertEqual(str(Sr5), "x*(%s + x*(%s + x))" % (f_1, f_1))
        self.assertEqual(str(Sr6), "x*x*(a + c + x*(b + d))")
        self.assertEqual(str(Sr7), "(x*(%s + x + y) + z*(%s + y))" % (f_2, f_2))
        self.assertEqual(str(Sr8), "(x*x*x*(y + z) + y*(a + b))")
        self.assertEqual(str(Sr9), "(x*x*(%s + x*(y + z)) + y*(a + b + c))" % f_2)
        self.assertEqual(str(Sr10), "x*x*y*(%s + z)" % f_2)
        self.assertEqual(str(Sr11), "x*x*y*y*(%s + z)" % f_2)
        self.assertEqual(str(Sr12), "(x*x*y*y*(%s + z) + z*(a + b + c))" % f_2)
        self.assertEqual(str(Sr13), "(%s/x + %s/y)" % (f_1, f_1))
        self.assertEqual(str(Sr14), "(-%s/x-%s/y)" % (f_1, f_1))
        self.assertEqual(str(Sr15), "%s/x" % format["float"](4))
        self.assertEqual(str(Sr16), "(%s + %s*x)/(y*z)" % (format["float"](0.5), f_2))
        self.assertEqual(str(Sr17), "x*y*(%s/a + z/b)" % f_2)
        self.assertEqual(str(Sr18), "(%s + x*(z + %s*y))/a" % (f_2, format["float"](3)))
        self.assertEqual(str(Sr19), "x*(z + (%s + y)/a)" % f_2)
        self.assertEqual(str(Sr20), "a*c*d*(x + y)/(b*z)")
        self.assertEqual(str(Sr21), "(x*(a + b + c + y + z) + y*(%s + a + b) + z*(%s + a + b))" % (f_2, f_2))
        self.assertEqual(str(Sr22), "%s" % format["float"](0))
        self.assertEqual(str(Sr23), "(x*y*z + z*z*z*(y*y*y*(x + z) + z*z))")

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
        self.assertEqual(S20.ops(), 6)
        self.assertEqual(Sr20.ops(), 6)
        self.assertEqual(S21.ops(), 21)
        self.assertEqual(Sr21.ops(), 13)
        self.assertEqual(S23.ops(), 21)
        self.assertEqual(Sr23.ops(), 12)

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestReduceOperations('testReduceOperations'))

