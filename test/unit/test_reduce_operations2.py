"Some simple functions for manipulating expressions symbolically"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-03-04 -- 2009-03-04"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"


import unittest
from ufl.common import dstr
from ffc.compiler.codegeneration.quadrature.reduce_operations import operation_count, expand_operations, reduce_operations
from ffc.compiler.codegeneration.quadrature.reduce_operations2 import *

import time

from ffc.compiler.format.ufcformat import Format
from ffc.common.constants import FFC_OPTIONS

class Tests(unittest.TestCase):

    def testSymbol(self):
        "Test simple symbol instance."

        s0 = Symbol("x", -1, BASIS)
        s1 = Symbol("y", 1, IP)
        s2 = Symbol("z", -2.1, GEO)
        s3 = Symbol("z", 3, GEO)
        s4 = Symbol("z", 0, IP)
#        print "\nTesting Symbols"
#        print "s0: '%s'" %s0
#        print "s1: '%s'" %s1
#        print "s2: '%s'" %s2
#        print "s3: '%s'" %s3
#        print "s4: '%s'" %s4
        self.assertEqual(s0.__repr__(), "Symbol('x', -1.0000, BASIS)")
        self.assertEqual(s1.__repr__(), "Symbol('y', 1.0000, IP)")
        self.assertEqual(s2.__repr__(), "Symbol('z', -2.1000, GEO)")
        self.assertEqual(s4.__repr__(), "Symbol('z', 0.0000, IP)")
        self.assertEqual(s2 == s3, True)
        self.assertEqual(s2 == s1, False)
        self.assertEqual(str(s0), ' - x')
        self.assertEqual(str(s1), 'y')
        self.assertEqual(str(s2), ' - 2.1*z')
        self.assertEqual(str(s3), '3*z')
        self.assertEqual(str(s4), '0')
        self.assertEqual(s0.ops(), 1)
        self.assertEqual(s1.ops(), 0)
        self.assertEqual(s2.ops(), 2)
        self.assertEqual(s3.ops(), 1)
        self.assertEqual(s4.ops(), 0)


    def testProduct(self):
        "Test simple product instance."

        x = Symbol("x", 2, GEO)
        y = Symbol("y", -2, IP)
        z = Symbol("z", 1, GEO)
        z0 = Symbol("z", 0, GEO)
        c = Symbol("", 4, CONST)

        p0 = Product([x, x], False)
        p1 = Product([x, y], False)
        p2 = Product([x, z, y], False)
        p3 = Product([y, None], False)
        p4 = Product([y, z0, x], False)
        p5 = Product([z], False)
        p6 = Product([x, c], False)

        s0 = Symbol(get_format()["cos"](p0), 2, p0.t)
        s0.base_expr = p0 # Set base expression
        s0.base_op = 1    # count cos() as one operation

#        print "\nTesting Products"
#        print "x:  '%s'" %x
#        print "y:  '%s'" %y
#        print "z:  '%s'" %z
#        print "z0: '%s'" %z0
#        print "\np0 = %s * %s = '%s'" % (x, x, p0)
#        print "\np1 = %s * %s = '%s'" % (x, y, p1)
#        print "\np2 = %s * %s * %s = '%s'" % (x, z, y, p2)
#        print "\np3 = %s * %s = '%s'" %(y, None, p3)
#        print "\np4 = %s * %s * %s = '%s'" %(y, z0, x, p4)
#        print "\np5 = %s = '%s'" % (z, p5)
#        print "\ns0 = 2 * cos(%s) = '%s'" % (p0, s0)
#        print "\np6 = %s * %s = '%s'" % (x, c, p6)

        self.assertEqual(p0.__repr__(), "Product([Symbol('', 4.0000, CONST), Symbol('x', 1.0000, GEO), Symbol('x', 1.0000, GEO)])")
        self.assertEqual(p0.t == GEO, True)
        self.assertEqual(p1.__repr__(), "Product([Symbol('', -4.0000, CONST), Symbol('y', 1.0000, IP), Symbol('x', 1.0000, GEO)])")
        self.assertEqual(p1.t == IP, True)
        self.assertEqual(p2.c < 0, True)
        self.assertEqual(p3.c == 0, True)
        self.assertEqual(p4.c == 0, True)
        self.assertEqual(p5.t == GEO, True)
        self.assertEqual(str(p0), '4*x*x')
        self.assertEqual(str(p1), ' - 4*y*x')
        self.assertEqual(str(p2), ' - 4*y*x*z')
        self.assertEqual(str(p3), '0')
        self.assertEqual(str(p4), '0')
        self.assertEqual(str(p5), 'z')
        self.assertEqual(str(p6), '8*x')
#        self.assertEqual(p0.ops(), 2)
#        self.assertEqual(p1.ops(), 3)
#        self.assertEqual(p2.ops(), 4)
#        self.assertEqual(p3.ops(), 0)
#        self.assertEqual(p4.ops(), 0)
#        self.assertEqual(p5.ops(), 0)
#        self.assertEqual(p6.ops(), 1)
#        self.assertEqual(str(s0), '2*std::cos(4*x*x)')
        self.assertEqual(s0.ops(), 4)


    def testSum(self):
        "Test simple sum instance."

        s0 = Symbol("x", 2, GEO)
        s1 = Symbol("y", 4.5, GEO)
        s2 = Symbol("z", -8, GEO)
        zp= Symbol("z", 4, GEO)
        c = Symbol("", 4, CONST)
        x0 = Symbol("x", 0, GEO)

        p0 = Product([s0, s1], False)
        p1 = Product([s0, s2], False)

        f0 = Fraction(s1, s0, False)
        f1 = Fraction(s0, s2, False)

        S0  = Sum([s0], False)
        xy  = Sum([s0, s1], False)
        xx  = Sum([s0, s0], False)
        zz  = Sum([s2, zp], False)
        xyz = Sum([s0, s1, zp], False)
        zzz  = Sum([s2, zp, zp], False)
        cx0z  = Sum([c, x0, s2], False)
        cNz  = Sum([c, None, s2], False)

        Sp0 = Sum([p0, p1], False)
        Sp1 = Sum([p1, p1], False)
        Sf0 = Sum([f0, f1], False)
        Sf1 = Sum([f1, f1], False)

#        print "\nTesting Sum"
#        print "%s  = '%s'" %(s0, S0)
#        print "Sx.recon(): '%s'" %(S0.recon())
#        print "%s + %s  = '%s'" %(s0, s1, xy)
#        print "%s + %s  = '%s'" %(s0, s0, xx)
#        print "%s + %s  = '%s'" %(s2, zp, zz)
#        print "%s + %s + %s = '%s'" %(s0, s1, zp, xyz)
#        print "%s + %s + %s = '%s'" %(s2, zp, zp, zzz)
#        print "%s + %s + %s = '%s'" %(c, x0, s0, cx0z)
#        print "%s + %s + %s = '%s'" %(c, None, s2, cNz)

#        print "\n%s * %s  = '%s'" %(s0, s1, p0)
#        print "%s * %s  = '%s'" %(s0, s2, p1)
#        print "%s + %s  = '%s'" %(p0, p1, Sp0)
#        print "%s + %s  = '%s'" %(p1, p1, Sp1)

#        print "\n%s / %s  = '%s'" %(s1, s0, f0)
#        print "%s / %s  = '%s'" %(s0, s2, f1)
#        print "%s + %s  = '%s'" %(f0, f1, Sf0)
#        print "%s + %s  = '%s'" %(f1, f1, Sf1)
        x = 2.3
        self.assertEqual(xx.__repr__(), "Sum([Symbol('x', 4.0000, GEO)])")
        self.assertEqual(str(xy), "(2*x + 4.5*y)")
        self.assertEqual(str(xx), "4*x")
        self.assertEqual(str(zz), " - 4*z")
        self.assertEqual(str(xyz), "(2*x + 4.5*y + 4*z)")
        self.assertEqual(str(zzz), "0")
        self.assertEqual(str(cx0z), "(4 - 8*z)")
        self.assertEqual(str(cNz), "(4 - 8*z)")
        self.assertEqual(str(Sp0), '(9*x*y - 16*x*z)')
        self.assertEqual(str(Sp1), ' - 32*x*z')
        self.assertEqual(str(Sf0), '(2.25*y/x - 0.25*x/z)')
        self.assertEqual(str(Sf1), ' - 0.5*x/z')
        self.assertEqual(eval(str(S0)), eval(str(s0)))
 
        self.assertEqual(xy.ops(), 3)
        self.assertEqual(xx.ops(), 1)
        self.assertEqual(zz.ops(), 2)
        self.assertEqual(xyz.ops(), 5)
        self.assertEqual(zzz.ops(), 0)
        self.assertEqual(cx0z.ops(), 2)
        self.assertEqual(cNz.ops(), 2)
        self.assertEqual(Sp0.ops(), 5)
        self.assertEqual(Sp1.ops(), 3)
        self.assertEqual(Sf0.ops(), 5)
        self.assertEqual(Sf1.ops(), 3)


    def testFraction(self):
        "Test simple fraction instance."

        x = Symbol("x", 2, GEO)
        y = Symbol("y", 4.5, GEO)
        z = Symbol("z", -8, GEO)
        zp= Symbol("z", 4, GEO)
        c0 = Symbol("", 4, CONST)
        c1 = Symbol("", 16, CONST)
        x0 = Symbol("x", 0, GEO)
        xp = Symbol("x", 1,  GEO)
        yp = Symbol("y", 1,  GEO)

        s0 = Sum([c0, x], False)
        s1 = Sum([x, x], False)
        p0 = Product([c0, x], False)
        p1 = Product([y, x], False)

        f0 = Fraction(y, x, False)
        f1 = Fraction(z, zp, False)
        f2 = Fraction(x0, y, False)
        f3 = Fraction(z, c0, False)
        f4 = Fraction(c0, x, False)
        f5 = Fraction(s0, p0, False)
        f6 = Fraction(p0, s0, False)
        f7 = Fraction(p0, f0, False)
        f8 = Fraction(xp, yp, False)
        f9 = Fraction(f8, p1, False)
        f10= Fraction(s1, c1, False)

#        print "\nTesting Fractions"
#        print "x:  '%s'" %x
#        print "y:  '%s'" %y
#        print "z:  '%s'" %z
#        print "zp: '%s'" %zp
#        print "c0:  '%s'" %c0
#        print "c1:  '%s'" %c1
#        print "x0: '%s'" %x0
#        print "f0 = frac(%s, %s) = '%s'" %(y, x, f0)
#        print "f1 = frac(%s, %s) = '%s'" %(z, zp, f1)
#        print "f2 = frac(%s, %s) = '%s'" %(x0, y, f2)
#        print "f3 = frac(%s, %s) = '%s'" %(z, c0, f3)
#        print "f4 = frac(%s, %s) = '%s'" %(c0, x, f4)
#        print "f5 = frac(%s, %s) = '%s'" %(s0, p0, f5)
#        print "f6 = frac(%s, %s) = '%s'" %(p0, s0, f6)
#        print "f7 = frac(%s, %s) = '%s'" %(p0, f0, f7)
#        print "f8 = frac(%s, %s) = '%s'" %(xp, yp, f8)
#        print "f9 = frac(%s, %s) = '%s'" %(f8, p1, f9)
#        print "f10= frac(%s, %s) = '%s'" %(s1, c1, f10)

        self.assertEqual(f0.__repr__(), "Fraction(Symbol('y', 2.2500, GEO), Symbol('x', 1.0000, GEO))")
        self.assertEqual(str(f0), "2.25*y/x")
        self.assertEqual(str(f1), " - 2")
        self.assertEqual(f1.t , CONST)
        self.assertEqual(str(f2), "0")
        self.assertEqual(str(f3), " - 2*z")
        self.assertEqual(str(f4), "2/x")
        self.assertEqual(str(f5), '0.125/x*(2*x + 4)')
        self.assertEqual(str(f6), '8*x/(2*x + 4)')
        self.assertEqual(str(f7), '3.55555555555556*x/(y/x)')
        self.assertEqual(str(f8), 'x/y')
        self.assertEqual(str(f9), '0.111111111111111*(x/y)/(x*y)')
        self.assertEqual(str(f10), '0.25*x')

        self.assertEqual(f0.ops(), 2)
        self.assertEqual(f1.ops(), 1)
        self.assertEqual(f2.ops(), 0)
        self.assertEqual(f3.ops(), 2)
        self.assertEqual(f4.ops(), 1)
        self.assertEqual(f5.ops(), 4)
        self.assertEqual(f6.ops(), 4)
        self.assertEqual(f7.ops(), 3)
        self.assertEqual(f8.ops(), 1)
        self.assertEqual(f9.ops(), 4)
        self.assertEqual(f10.ops(), 1)

    def testMixedSymbols(self):

        s0 = Symbol("x", -1, GEO)
        s1 = Symbol("y", 1, GEO)
        s2 = Symbol("z", 3, GEO)
        s3= Symbol("z", -4, GEO)

#        print "\nTesting Mixed Symbols"
#        print "s0: '%s'" %s0
#        print "s1: '%s'" %s1
#        print "s2: '%s'" %s2
#        print "s3: '%s'" %s3

        x = 2.2
        y = -82.2
        z = 10.2

        e0 = Product([s0, Sum([s0, s3], False)], False)
        e1 = Fraction(Product([s0, Sum([s0, s3], False)], False), s3, False)
        e2 = Fraction(Product([s0, Sum([s0, s3], False)], False), Sum([s0, s3], False), False)
        e3 = Product([s0, Sum([s0, s2], False), Sum([s0, s1], False)], False)
        e4 = Product([e0, e0, e3, s0], False)
        e5 = Sum([e0, e0, e3, s0], False)
        e6 = Sum([s1, Sum([s0, s3], False), Sum([s0, s3], False)], False)

#        print "e0: '%s'" %e0
#        print "e1: '%s'" %e1
#        print "e2: '%s'" %e2
#        print "e3: '%s'" %e3
#        print "e4: '%s'" %e4
#        print "e5: '%s'" %e5
#        print "e6: '%s'" %e6
        self.assertEqual(str(e0), 'x*(x + 4*z)')
        self.assertEqual(str(e1), ' - 0.25/z*x*(x + 4*z)')
        self.assertEqual(str(e2), ' - x/(x + 4*z)*(x + 4*z)')
        self.assertEqual(str(e3), ' - x*(y - x)*(3*z - x)')
        self.assertEqual(str(e4), 'x*x*(x + 4*z)*x*(x + 4*z)*x*(y - x)*(3*z - x)')
        self.assertEqual(str(e5), '(2*x*(x + 4*z) - x - x*(y - x)*(3*z - x))')
        self.assertAlmostEqual(eval(str(e0)), eval(str(s0))*(eval(str(s0)) + eval(str(s3))))
        self.assertAlmostEqual(eval(str(e1)), eval(str(s0))*(eval(str(s0)) + eval(str(s3)))/eval(str(s3)))
        self.assertAlmostEqual(eval(str(e2)), eval(str(s0))*(eval(str(s0)) + eval(str(s3)))/(eval(str(s0)) + eval(str(s3))))
        self.assertAlmostEqual(eval(str(e3)), eval(str(s0))*(eval(str(s0)) + eval(str(s2)))*(eval(str(s0)) + eval(str(s1))))
        self.assertAlmostEqual(eval(str(e4)), eval(str(e0))*eval(str(e0))*eval(str(e3))*eval(str(s0)))
        self.assertAlmostEqual(eval(str(e5)), eval(str(e0))+eval(str(e0))+eval(str(e3))+eval(str(s0)))
        self.assertAlmostEqual(eval(str(e6)), eval(str(s1))+eval(str(s0)) + eval(str(s3))+eval(str(s0)) + eval(str(s3)))
        self.assertEqual(e0.ops(), 3)
        self.assertEqual(e1.ops(), 6)
        self.assertEqual(e2.ops(), 7)
        self.assertEqual(e3.ops(), 6)
        self.assertEqual(e4.ops(), 14)
        self.assertEqual(e5.ops(), 11)
        self.assertEqual(e6.ops(), 4)

    def testRemoveNested(self):

        s0 = Symbol("x", -1, GEO)
        s1 = Symbol("y", 1, GEO)
        s2 = Symbol("z", 5, IP)
        s3 = Symbol("z", -4, GEO)
        s4 = Symbol("w", 0, GEO)

        x = 2.2
        y = -0.2
        z = 1.1

#        print "\nTesting Nested Symbols"
#        print "s0: '%s'" %s0
#        print "s1: '%s'" %s1
#        print "s2: '%s'" %s2
#        print "s3: '%s'" %s3
#        print "s4: '%s'" %s4

        P0 = Product([s2, s1], False)
        P1 = Product([P0, s0], False)
        P2 = Product([P1, s1, P0], False)
        P3 = Product([P1, P2], False)

        S0 = Sum([s2, s1], False)
        S1 = Sum([S0, s0], False)
        S2 = Sum([S1, s1, S0], False)
        S3 = Sum([S1, S2], False)
        S4 = Sum([s0], False)
        S4.remove_nested()

        F0 = Fraction(s2, s1, False)
        F1 = Fraction(F0, s0, False)
        F2 = Fraction(F1, F0, False)
        F3 = Fraction(F1, F2, False)

        # Special tests for Fraction
        F4 = Fraction(P0, F0, False)
        F4n = F4.remove_nested()
        F5 = Fraction(Fraction(s0, P0, False), P0, False)
        F5n = F5.remove_nested()
        F6 = Fraction( Fraction( Fraction(s1, s0, False), Fraction(s1, s2, False), False), Fraction( Fraction(s2, s0, False), Fraction(s1, s0, False), False), False )
        F6n = F6.remove_nested()
#        print "\nF4: '%s'" %F4
#        print "F4n: '%s'" %F4n
#        print "\nF4: '%s'" %F4.__repr__()
#        print "F4n: '%s'" %F4n.__repr__()
#        print "\nF5: '%s'" %F5
#        print "F5n: '%s'" %F5n
#        print "\nF5: '%s'" %F5.__repr__()
#        print "F5n: '%s'" %F5n.__repr__()
#        print "\nF6: '%s'" %F6
#        print "F6n: '%s'" %F6n
#        print "\nF6: '%s'" %F6.__repr__()
#        print "F6n: '%s'" %F6n.__repr__()

        e0 = Product([P3, F2], False)
        e1 = Product([S3, P2], False)
        e2 = Product([F3, S1], False)

        e3 = Sum([P3, F2], False)
        e4 = Sum([S3, P2], False)
        e5 = Sum([F3, S1], False)

        e6 = Fraction(P3, F2, False)
        e7 = Fraction(S3, P2, False)
        e8 = Fraction(F3, S1, False)

        ex0 = e0.remove_nested()
        ex1 = e1.remove_nested()
        ex2 = e2.remove_nested()
        ex3 = e3.remove_nested()
        ex4 = e4.remove_nested()
        ex5 = e5.remove_nested()
        ex6 = e6.remove_nested()
        ex7 = e7.remove_nested()
        ex8 = e8.remove_nested()
#        print "\ne0: '%s'" %e0
#        print "ex0: '%s'" %ex0
#        print "\ne1: '%s'" %e1
#        print "ex1: '%s'" %ex1
#        print "\ne2: '%s'" %e2
#        print "ex2: '%s'" %ex2
#        print "\ne3: '%s'" %e3
#        print "ex3: '%s'" %ex3
#        print "\ne4: '%s'" %e4
#        print "ex4: '%s'" %ex4
#        print "\ne5: '%s'" %e5
#        print "ex5: '%s'" %ex5
#        print "\ne6: '%s'" %e6
#        print "ex6: '%s'" %ex6
#        print "\ne7: '%s'" %e7
#        print "ex7: '%s'" %ex7
#        print "\ne8: '%s'" %e8
#        print "ex8: '%s'" %ex8

        self.assertAlmostEqual(eval(str(e0)), eval(str(ex0)))
        self.assertAlmostEqual(eval(str(e1)), eval(str(ex1)))
        self.assertAlmostEqual(eval(str(e2)), eval(str(ex2)))
        self.assertAlmostEqual(eval(str(e3)), eval(str(ex3)))
        self.assertAlmostEqual(eval(str(e4)), eval(str(ex4)))
        self.assertAlmostEqual(eval(str(e5)), eval(str(ex5)))
        self.assertAlmostEqual(eval(str(e6)), eval(str(ex6)))
        self.assertAlmostEqual(eval(str(e7)), eval(str(ex7)))
        self.assertAlmostEqual(eval(str(e8)), eval(str(ex8)))
        self.assertAlmostEqual(eval(str(F4)), eval(str(F4n)))
        self.assertAlmostEqual(eval(str(F5)), eval(str(F5n)))
        self.assertAlmostEqual(eval(str(F6)), eval(str(F6n)))

        self.assertEqual(e0.ops(), 15)
        self.assertEqual(ex0.ops(), 15)
        self.assertEqual(e1.ops(), 19)
        self.assertEqual(ex1.ops(), 13)
        self.assertEqual(e2.ops(), 12)
        self.assertEqual(ex2.ops(), 12)
        self.assertEqual(e3.ops(), 14)
        self.assertEqual(ex3.ops(), 14)
        self.assertEqual(e4.ops(), 18)
        self.assertEqual(ex4.ops(), 12)
        self.assertEqual(e5.ops(), 12)
        self.assertEqual(ex5.ops(), 12)
        self.assertEqual(e6.ops(), 15)
        self.assertEqual(ex6.ops(), 15)
        self.assertEqual(e7.ops(), 19)
        self.assertEqual(ex7.ops(), 13)
        self.assertEqual(e8.ops(), 12)
        self.assertEqual(ex8.ops(), 12)
        self.assertEqual(F4.ops(), 3)
        self.assertEqual(F4n.ops(), 3)
        self.assertEqual(F5.ops(), 6)
        self.assertEqual(F5n.ops(), 6)
        self.assertEqual(F6.ops(), 8)
        self.assertEqual(F6n.ops(), 8)

    def testExpandOperations(self):

        s0 = Symbol("x", -1, GEO)
        s1 = Symbol("y", 1, GEO)
        s2 = Symbol("z", 5, IP)
        s3 = Symbol("z", -4, GEO)
        s4 = Symbol("w", 0, GEO)
        c0 = Symbol("", 2, CONST)
        c1 = Symbol("", 1, CONST)

        x = 2.2
        y = -0.2
        z = 1.1

#        print "\nTesting Expanded"
#        print "s0: '%s'" %s0
#        print "s1: '%s'" %s1
#        print "s2: '%s'" %s2
#        print "s3: '%s'" %s3
#        print "s4: '%s'" %s4

        P0 = Product([s2, s1], False)
        P1 = Product([P0, s0], False)
        P2 = Product([P1, s1, P0], False)
        P3 = Product([P1, P2], False)

        S0 = Sum([s2, s1], False)
        S1 = Sum([S0, s0], False)
        S2 = Sum([S1, s1, S0], False)
        S3 = Sum([S1, S2], False)

        F0 = Fraction(s2, s1, False)
        F1 = Fraction(F0, s0, False)
        F2 = Fraction(F1, F0, False)
        F3 = Fraction(F1, F2, False)
        PF = Product([F0, F1], False)

        e0 = Product([P3, F2], False)
        e1 = Product([S3, P2], False)
        e2 = Product([F3, S1], False)
        e2p = Product([s1, S1], False)
        e2p.c = -1

        e3 = Sum([P3, F2], False)
        e4 = Sum([S3, P2], False)
        e5 = Sum([F3, S1], False)
        e5p = Sum([P0, e2p], False)

        e6 = Fraction(P3, F2, False)
        e7 = Fraction(S3, P2, False)
        e8 = Fraction(F3, S1, False)
        e9 = Fraction(S0, s0, False)
        e9 = Fraction(S0, s0, False)

        P00 = Product([Symbol("Jinv_00", 1, GEO)]*2, False)
        P01 = Product([Symbol("Jinv_01", 1, GEO)]*2, False)
        P20 = Product([Symbol("Jinv_00", 1, GEO), Symbol("Jinv_20", 2, GEO)])
        P21 = Product([Symbol("Jinv_01", 1, GEO), Symbol("Jinv_21", 2, GEO)])
        PS0 = Product([Symbol("Jinv_22", 1, GEO), Sum([P00, P01], False)])
        PS1 = Product([Symbol("Jinv_02", -1, GEO), Sum([P20, P21], False)])
        SP = Sum([PS0, PS1], False)

        E0 = Fraction(Sum([Product([c0])]), c1)

        ex0 = e0.expand()
        ex1 = e1.expand()
        ex2 = e2.expand()
        ex2p = e2p.expand()
        ex3 = e3.expand()
        ex4 = e4.expand()
        ex5 = e5.expand()
        ex5p = e5p.expand()
        ex6 = e6.expand()
        ex7 = e7.expand()
        ex8 = e8.expand()
        ex9 = e9.expand()
        PF.expand()
        SPx = SP.expand()
        Ex0 = E0.expand()
#        print "\ne0: '%s'" %e0
#        print "ex0: '%s'" %ex0
#        print "\ne1: '%s'" %e1
#        print "ex1: '%s'" %ex1
#        print "\ne2: '%s'" %e2
#        print "ex2: '%s'" %ex2
#        print "\ne2p: '%s'" %e2p
#        print "ex2p: '%s'" %ex2p
#        print "\ne3: '%s'" %e3
#        print "ex3: '%s'" %ex3
#        print "\ne4: '%s'" %e4
#        print "ex4: '%s'" %ex4
#        print "\ne5: '%s'" %e5
#        print "ex5: '%s'" %ex5
#        print "\ne5p: '%s'" %e5p
#        print "ex5p: '%s'" %ex5p
#        print "\ne6: '%s'" %e6
#        print "ex6: '%s'" %ex6
#        print "\ne7: '%s'" %e7
#        print "ex7: '%s'" %ex7
#        print "\ne8: '%s'" %e8
#        print "ex8: '%s'" %ex8
#        print "\ne9: '%s'" %e9
#        print "ex9: '%s'" %ex9
#        print "\nSP: '%s'" %SP
#        print "SPx: '%s'" %SPx
#        print "\nE0: '%s'" %E0
#        print "Ex0: '%s'" %Ex0

        Jinv_00, Jinv_01, Jinv_10, Jinv_02, Jinv_20, Jinv_22, Jinv_21, W1, det = [1,2,3,4,5,6,7,8,9]

        self.assertAlmostEqual(eval(str(e0)), eval(str(e0.remove_nested())))
        self.assertAlmostEqual(eval(str(e0)), eval(str(ex0)))
        self.assertAlmostEqual(eval(str(e1)), eval(str(ex1)))
        self.assertAlmostEqual(eval(str(e2)), eval(str(ex2)))
        self.assertAlmostEqual(eval(str(e3)), eval(str(ex3)))
        self.assertAlmostEqual(eval(str(e4)), eval(str(ex4)))
        self.assertAlmostEqual(eval(str(e5)), eval(str(ex5)))
        self.assertAlmostEqual(eval(str(e6)), eval(str(ex6)))
        self.assertAlmostEqual(eval(str(e7)), eval(str(ex7)))
        self.assertAlmostEqual(eval(str(e8)), eval(str(ex8)))
        self.assertAlmostEqual(eval(str(e9)), eval(str(ex9)))
        self.assertAlmostEqual(eval(str(SP)), eval(str(SPx)))
        self.assertAlmostEqual(eval(str(E0)), eval(str(Ex0)))

        self.assertEqual(e0.ops(), 15)
        self.assertEqual(ex0.ops(), 9)
        self.assertEqual(e1.ops(), 19)
        self.assertEqual(ex1.ops(), 23)
        self.assertEqual(e2.ops(), 12)
        self.assertEqual(ex2.ops(), 9)
        self.assertEqual(e3.ops(), 14)
        self.assertEqual(ex3.ops(), 11)
        self.assertEqual(e4.ops(), 18)
        self.assertEqual(ex4.ops(), 12)
        self.assertEqual(e5.ops(), 12)
        self.assertEqual(ex5.ops(), 6)
        self.assertEqual(e6.ops(), 15)
        self.assertEqual(ex6.ops(), 11)
        self.assertEqual(e7.ops(), 19)
        self.assertEqual(ex7.ops(), 17)
        self.assertEqual(e8.ops(), 12)
        self.assertEqual(ex8.ops(), 8)
        self.assertEqual(e9.ops(), 5)
        self.assertEqual(ex9.ops(), 5)
        self.assertEqual(Ex0.ops(), 0)

    def testReduceVarType(self):

        B0 = Symbol("B0", 1, BASIS)
        B1 = Symbol("B1", -4, BASIS)
        I0 = Symbol("I0", 5, IP)
        I1 = Symbol("I0", 1, IP)
        G0 = Symbol("G0", 3, GEO)
        C0 = Symbol("C0", 2, CONST)
        one = Symbol("", 1, CONST)
        five = Symbol("", 5, CONST)
        mfour = Symbol("", -4, CONST)

        p0 = Product([B0,I0], False)
        p1 = Product([B0,B1], False)

        s0 = Sum([B0,I0], False)
        s1 = Sum([p0,p1], False)
        s2 = Sum([B0,B1], False)
        s3 = Sum([B0,p0], False)

        f0 = Fraction(B0,I0, False).expand()
        f1 = Fraction(p1,I0, False).expand()
        f2 = Fraction(G0,s2, False).expand()
        f3 = Fraction(G0,s3, False).expand()

        r0 = B0.reduce_vartype(BASIS)
        r1 = B0.reduce_vartype(CONST)

        rp0 = p0.reduce_vartype(BASIS)
        rp1 = p0.reduce_vartype(IP)
        rp2 = p1.reduce_vartype(BASIS)
        rp3 = p1.reduce_vartype(GEO)

        rs0 = s0.reduce_vartype(BASIS)
        rs1 = s0.reduce_vartype(IP)
        rs2 = s1.reduce_vartype(BASIS)

        rf0 = f0.reduce_vartype(BASIS)
        rf1 = f1.reduce_vartype(BASIS)
        rf2 = f0.reduce_vartype(IP)
        rf3 = f2.reduce_vartype(BASIS)
        rf4 = f3.reduce_vartype(BASIS)


#        print "%s, red(BASIS): ('%s', '%s')" %(B0, r0[0], r0[1])
#        print "%s, red(CONST): ('%s', '%s')" %(B0, r1[0], r1[1])

#        print "%s, red(BASIS): ('%s', '%s')" %(p0, rp0[0], rp0[1])
#        print "%s, red(IP):    ('%s', '%s')" %(p0, rp1[0], rp1[1])
#        print "%s, red(BASIS): ('%s', '%s')" %(p1, rp2[0], rp2[1])
#        print "%s, red(CONST): ('%s', '%s')" %(p1, rp3[0], rp3[1])

#        print "%s, red(BASIS): ('%s', '%s')" %(s0, rs0[0], rs0[1])
#        print "%s, red(IP):    ('%s', '%s')" %(s0, rs1[0], rs1[1])
#        print "%s, red(BASIS): '%s', '%s'" %(s1, rs2[0], rs2[1])

#        print "%s, red(BASIS): ('%s', '%s')" %(f0, rf0[0], rf0[1])
#        print "%s, red(BASIS): ('%s', '%s')" %(f1, rf1[0], rf1[1])
#        print "%s, red(IP): ('%s', '%s')" %(f0, rf2[0], rf2[1])
#        print "%s, red(BASIS): ('%s', '%s')" %(f2, rf3[0], rf3[1])
#        print "%s, red(BASIS): ('%s', '%s')" %(f3, rf4[0], rf4[1])

        self.assertEqual((B0, one), r0)
        self.assertEqual(([], B0), r1)

        self.assertEqual((B0, I0), rp0)
        self.assertEqual((I0, B0),  rp1)
        self.assertEqual((p1, one), rp2)
        self.assertEqual(([], p1),  rp3)

        self.assertEqual(((), I0), rs0[0])
        self.assertEqual((B0, one), rs0[1])
        self.assertEqual((I1, five), rs1[0])
        self.assertEqual(((), B0), rs1[1])
        self.assertEqual((Product([Symbol("B0", 1, BASIS), Symbol("B1", 1, BASIS)]), mfour), rs2[0])
        self.assertEqual((B0, I0), rs2[1])

        self.assertEqual((B0, Fraction(Symbol("", 0.2, CONST), I0)), rf0)
        self.assertEqual((Product([Symbol("B0", 1, BASIS), Symbol("B1", 1, BASIS)]), Fraction(Symbol("", -0.8, CONST), Symbol("I0", 1, IP))), rf1)
        self.assertEqual( ( Fraction(Symbol("", 1.0, CONST), Symbol("I0", 1, IP)), Symbol("B0", 0.2, BASIS)), rf2)
        self.assertEqual(([], f2), rf3)
        self.assertEqual( ( Fraction(Symbol("", 1.0, CONST), Symbol("B0", 1, BASIS)),
          Fraction( Symbol("G0", 3.0, GEO), Sum([I0, one]))), rf4)


    def testReduceOperations(self):

        s0 = Symbol("x", -1, GEO)
        s1 = Symbol("y", 1, GEO)
        s2 = Symbol("z", 5, IP)
        s3 = Symbol("z", -4, GEO)
        s4 = Symbol("w", 0, GEO)

        x = 2.2
        y = -0.2
        z = 1.1

#        print "\nTesting ReduceOperations"
#        print "s0: '%s'" %s0
#        print "s1: '%s'" %s1
#        print "s2: '%s'" %s2
#        print "s3: '%s'" %s3
#        print "s4: '%s'" %s4

        P0 = Product([s2, s1], False)
        P1 = Product([P0, s0], False)
        P2 = Product([P1, s1, P0], False)
        P3 = Product([P1, P2], False)

        S0 = Sum([s2, s1], False)
        S1 = Sum([S0, s0], False)
        S2 = Sum([S1, s1, S0], False)
        S3 = Sum([S1, S2], False)

        F0 = Fraction(s2, s1, False)
        F1 = Fraction(F0, s0, False)
        F2 = Fraction(F1, F0, False)
        F3 = Fraction(F1, F2, False)

        e0 = Product([P3, F2], False)
        e1 = Product([S3, P2], False)
        e2 = Product([F3, S1], False)

        e3 = Sum([P3, F2], False)
        e4 = Sum([S3, P2], False)
        e5 = Sum([F3, S1], False)

        e6 = Fraction(P3, F2, False)
        e7 = Fraction(S3, P2, False)
        e8 = Fraction(F3, S1, False)
        e9 = Sum([Product([s0, s0, Sum([s1, Symbol("", 1, CONST)])]), Fraction(Sum([s1, Symbol("", 1, CONST)]), Product([s0,s0]))])

        ex0 = e0.expand()
        ex1 = e1.expand()
        ex2 = e2.expand()
        ex3 = e3.expand()
        ex4 = e4.expand()
        ex5 = e5.expand()
        ex6 = e6.expand()
        ex7 = e7.expand()
        ex8 = e8.expand()
        ex9 = e9.expand()

        er0 = ex0.reduce_ops()
        er1 = ex1.reduce_ops()
        er2 = ex2.reduce_ops()
        er3 = ex3.reduce_ops()
        er4 = ex4.reduce_ops()
        er5 = ex5.reduce_ops()
        er6 = ex6.reduce_ops()
#        er7 = ex7.reduce_ops()
        er8 = ex8.reduce_ops()
        er9 = ex9.reduce_ops()

#        print "\nex0: '%s'" %ex0
#        print "er0: '%s'" %er0
#        print "\nex1: '%s'" %ex1
#        print "er1: '%s'" %er1
#        print "\nex2: '%s'" %ex2
#        print "er2: '%s'" %er2
#        print "\nex3: '%s'" %ex3
#        print "er3: '%s'" %er3
#        print "\nex4: '%s'" %ex4
#        print "er4: '%s'" %er4
#        print "\nex5: '%s'" %ex5
#        print "er5: '%s'" %er5
#        print "\nex6: '%s'" %ex6
#        print "er6: '%s'" %er6
#        print "\nex7: '%s'" %ex7
#        print "er7: '%s'" %er7
#        print "\nex8: '%s'" %ex8
#        print "er8: '%s'" %er8
#        print "\nex9: '%s'" %ex9
#        print "er9: '%s'" %er9

#        self.assertAlmostEqual(eval(str(e0)), eval(str(e0.remove_nested())))
#        self.assertAlmostEqual(eval(str(ex0)), eval(str(er0)))
#        self.assertAlmostEqual(eval(str(ex1)), eval(str(er1)))
#        self.assertAlmostEqual(eval(str(ex2)), eval(str(er2)))
#        self.assertAlmostEqual(eval(str(ex3)), eval(str(er3)))
#        self.assertAlmostEqual(eval(str(ex4)), eval(str(er4)))
#        self.assertAlmostEqual(eval(str(ex5)), eval(str(er5)))
#        self.assertAlmostEqual(eval(str(ex6)), eval(str(er6)))
#        self.assertAlmostEqual(eval(str(ex7)), eval(str(er7)))
#        self.assertAlmostEqual(eval(str(ex8)), eval(str(er8)))
#        self.assertAlmostEqual(eval(str(ex9)), eval(str(er9)))

#        self.assertEqual(ex0.ops(), 9)
#        self.assertEqual(er0.ops(), 9)
#        self.assertEqual(ex1.ops(), 23)
#        self.assertEqual(er1.ops(), 11)
#        self.assertEqual(ex2.ops(), 9)
#        self.assertEqual(er2.ops(), 7)
#        self.assertEqual(ex3.ops(), 11)
#        self.assertEqual(er3.ops(), 11)
#        self.assertEqual(ex4.ops(), 12)
#        self.assertEqual(er4.ops(), 12)
#        self.assertEqual(ex5.ops(), 6)
#        self.assertEqual(er5.ops(), 6)
#        self.assertEqual(ex6.ops(), 11)
#        self.assertEqual(er6.ops(), 11)
#        self.assertEqual(ex7.ops(), 17)
#        self.assertEqual(er7.ops(), 17)
#        self.assertEqual(ex8.ops(), 8)
#        self.assertEqual(er8.ops(), 6)
#        self.assertEqual(ex9.ops(), 10)
#        self.assertEqual(er9.ops(), 8)


    def testDGElastoDyn(self):
        expr = Product([
                       Sum([
                            Symbol("F0", 1, IP),
                            Symbol("F1",-1,IP)
                          ]),
                       Fraction(
                                 Symbol("w4", 1, GEO),
                                 Symbol("w3", 1, GEO)
                                ),
                       Fraction(
                                 Product([
                                          Symbol("w2", 1, GEO),
                                          Symbol("w5", 1, GEO)
                                         ]),
                                 Symbol("w6", 1, GEO)
                                )
                      ])

        print
        start = time.time()
        expr_rem = expr.remove_nested()
        print "DGElastoDyn: time, remove_nested():                  ", time.time() - start

        start = time.time()
        expr_exp = expr.expand()
        print "DGElastoDyn: time, expand() (incl. remove_nested()): ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "DGElastoDyn: time, reduce_ops():                     ", time.time() - start

        print "expr.ops():     ", expr.ops()
        print "expr_exp.ops(): ", expr_exp.ops()
        print "expr_red.ops(): ", expr_red.ops()

#        print "expr:\n", expr
#        print "exp:\n", expr_exp
#        print "red:\n", expr_red

        F0, F1, w2, w3, w4, w5, w6 = (3.12, -8.1, -45.3, 17.5, 2.2, 5.3, 9.145)
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_rem)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))

    def testReduceGIP(self):

        expr = Sum([
                    Product([
                              Symbol("F17", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP), Symbol("F18", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G0", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F13", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G1", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F13", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F20", 1, IP),
                              Symbol("F3", 1, IP), Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F18", 1, IP), Symbol("F20", 1, IP), Symbol("F3", 1, IP),
                              Symbol("F8", 1, IP), Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G3", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP), Symbol("F20", 1, IP),
                              Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F20", 1, IP), Symbol("F8", 1, IP), Symbol("F8", 1, IP),
                              Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ]),
                    Product([
                              Symbol("F20", 1, IP), Symbol("F8", 1, IP), Symbol("F8", 1, IP), Symbol("F9", 1, IP),
                              Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G5", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G6", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F10", 1, IP), Symbol("F20", 1, IP), Symbol("F3", 1, IP),
                              Symbol("F8", 1, IP), Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G1", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F10", 1, IP), Symbol("F20", 1, IP), Symbol("F8", 1, IP),
                              Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G7", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP), Symbol("F19", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F17", 1, IP), Symbol("F19", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F13", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G7", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F17", 1, IP), Symbol("F19", 1, IP), Symbol("F20", 1, IP),
                              Symbol("F3", 1, IP), Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G8", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP), Symbol("F18", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G5", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F17", 1, IP), Symbol("F19", 1, IP), Symbol("F20", 1, IP),
                              Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G9", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP), Symbol("F20", 1, IP),
                              Symbol("F3", 1, IP), Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F18", 1, IP), Symbol("F20", 1, IP), Symbol("F8", 1, IP),
                              Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G6", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F20", 1, IP),
                              Symbol("F3", 1, IP), Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G8", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F13", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F20", 1, IP),
                              Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G9", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F20", 1, IP), Symbol("F3", 1, IP),
                              Symbol("F8", 1, IP), Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G3", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("F8", 1, IP),
                              Symbol("F8", 1, IP), Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F19", 1, IP), Symbol("F20", 1, IP), Symbol("F8", 1, IP),
                              Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F17", 1, IP), Symbol("F19", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G9", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F12", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G0", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F17", 1, IP), Symbol("F19", 1, IP), Symbol("F19", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G7", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F17", 1, IP), Symbol("F18", 1, IP), Symbol("F19", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("F8", 1, IP), Symbol("F8", 1, IP),
                              Symbol("F9", 1, IP), Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G0", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F20", 1, IP), Symbol("F8", 1, IP),
                              Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G6", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F17", 1, IP), Symbol("F19", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G8", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F19", 1, IP), Symbol("F20", 1, IP), Symbol("F3", 1, IP),
                              Symbol("F8", 1, IP), Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ]),
                    Product([
                              Symbol("F10", 1, IP), Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F20", 1, IP),
                              Symbol("F8", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F13", 1, IP), Symbol("F20", 1, IP), Symbol("F8", 1, IP),
                              Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G4", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F12", 1, IP),
                              Symbol("F20", 1, IP), Symbol("W9", 1, IP), Symbol("G5", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F20", 1, IP), Symbol("F3", 1, IP),
                              Symbol("F8", 1, IP), Symbol("F9", 1, IP), Symbol("W9", 1, IP), Symbol("G3", 1, GEO)
                            ]),
                    Product([
                              Symbol("F17", 1, IP), Symbol("F17", 1, IP), Symbol("F19", 1, IP), Symbol("F19", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G1", 1, GEO)
                            ]),
                    Product([
                              Symbol("F11", 1, IP), Symbol("F12", 1, IP), Symbol("F17", 1, IP), Symbol("F19", 1, IP),
                              Symbol("F20", 1, IP), Symbol("F3", 1, IP), Symbol("W9", 1, IP), Symbol("G2", 1, GEO)
                            ])
                   ])

        print
        start = time.time()
        expr_rem = expr.remove_nested()
        print "ReduceGIP: time, remove_nested():                  ", time.time() - start

        start = time.time()
        expr_exp = expr.expand()
        print "ReduceGIP: time, expand() (incl. remove_nested()): ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "ReduceGIP: time, reduce_ops():                     ", time.time() - start

        print "expr.ops():     ", expr.ops()
        print "expr_exp.ops(): ", expr_exp.ops()
        print "expr_red.ops(): ", expr_red.ops()

#        print "expr: ", expr
#        print "exp:  ", expr_exp
#        print "red:  ", expr_red

        W9 = 9
        F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20 = [0.123 * i for i in range(1,21)]
        G0, G1, G2, G3, G4, G5, G6, G7, G8, G9 = [2.64 + 1.0/i for i in range(20, 30)]

        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_rem)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(expr.ops() > expr_red.ops(), True)


    def testPoisson(self):

        poisson = """((Jinv_00*FE0_D10_ip_j + Jinv_10*FE0_D01_ip_j)*(Jinv_00*FE0_D10_ip_k + Jinv_10*FE0_D01_ip_k) + (Jinv_01*FE0_D10_ip_j + Jinv_11*FE0_D01_ip_j)*(Jinv_01*FE0_D10_ip_k + Jinv_11*FE0_D01_ip_k))*W4_ip*det"""

        expr = Product([
                     Sum([
                          Product([
                                   Sum([
                                        Product([Symbol("Jinv_00", 1, GEO), Symbol("FE0_D10_ip_j", 1, BASIS)])
                                        ,
                                        Product([Symbol("Jinv_10", 1, GEO), Symbol("FE0_D01_ip_j", 1, BASIS)])
                                       ]),
                                   Sum([
                                        Product([Symbol("Jinv_00", 1, GEO), Symbol("FE0_D10_ip_k", 1, BASIS)])
                                        ,
                                        Product([Symbol("Jinv_10", 1, GEO), Symbol("FE0_D01_ip_k", 1, BASIS)])
                                       ])
                                  ])
                          ,
                          Product([
                                   Sum([
                                        Product([Symbol("Jinv_01", 1, GEO), Symbol("FE0_D10_ip_j", 1, BASIS)])
                                        ,
                                        Product([Symbol("Jinv_11", 1, GEO), Symbol("FE0_D01_ip_j", 1, BASIS)])
                                       ]),
                                   Sum([
                                        Product([Symbol("Jinv_01", 1, GEO), Symbol("FE0_D10_ip_k", 1, BASIS)])
                                        ,
                                        Product([Symbol("Jinv_11", 1, GEO), Symbol("FE0_D01_ip_k", 1, BASIS)])
                                       ])
                                  ])
                         ])
                     ,
                     Symbol("W4_ip", 1, IP)
                     ,
                     Symbol("det", 1, GEO)
                    ])

        print
        start = time.time()
        expr_rem = expr.remove_nested()
        print "Poisson: time, remove_nested():                  ", time.time() - start

        start = time.time()
        expr_exp = expr.expand()
        print "Poisson: time, expand() (incl. remove_nested()): ", time.time() - start

        start = time.time()
        poisson_exp = expand_operations(poisson, get_format())
        print "Poisson: time, old expand():                     ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "Poisson: time, reduce_ops():                     ", time.time() - start

        start = time.time()
        poisson_red = reduce_operations(poisson, get_format())
        print "Poisson: time, old reduce():                     ", time.time() - start

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
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_rem)))
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
                                             Product([Symbol("Jinv_00", 1, GEO), Symbol("FE0_C0_D10_ip_j", 1, BASIS)])
                                             ,
                                             Product([Symbol("Jinv_10", 1, GEO), Symbol("FE0_C0_D01_ip_j", 1, BASIS)])
                                            ])
                                        ,
                                        Symbol("", 2, CONST)
                                        ,
                                        Sum([
                                             Product([Symbol("Jinv_00", 1, GEO), Symbol("FE0_C0_D10_ip_k", 1, BASIS)])
                                             ,
                                             Product([Symbol("Jinv_10", 1, GEO), Symbol("FE0_C0_D01_ip_k", 1, BASIS)])
                                            ])
                                        ,
                                        Symbol("", 2, CONST)
                                        ])
                               ,
                               Product([
                                        Sum([
                                             Sum([
                                                  Product([Symbol("Jinv_00", 1, GEO), Symbol("FE0_C1_D10_ip_j", 1, BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_10", 1, GEO), Symbol("FE0_C1_D01_ip_j", 1, BASIS)])
                                                 ])
                                             ,
                                             Sum([
                                                  Product([Symbol("Jinv_01", 1, GEO), Symbol("FE0_C0_D10_ip_j", 1, BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_11", 1, GEO), Symbol("FE0_C0_D01_ip_j", 1, BASIS)])
                                                 ])
                                            ])
                                        ,
                                        Sum([
                                             Sum([
                                                  Product([Symbol("Jinv_00", 1, GEO), Symbol("FE0_C1_D10_ip_k", 1, BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_10", 1, GEO), Symbol("FE0_C1_D01_ip_k", 1, BASIS)])
                                                 ])
                                             ,
                                             Sum([
                                                  Product([Symbol("Jinv_01", 1, GEO), Symbol("FE0_C0_D10_ip_k", 1, BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_11", 1, GEO), Symbol("FE0_C0_D01_ip_k", 1, BASIS)])
                                                 ])
                                            ])
                                       ])
                              ])
                          ,
                          Sum([
                               Product([
                                        Sum([
                                             Product([Symbol("Jinv_01", 1, GEO), Symbol("FE0_C1_D10_ip_j", 1, BASIS)])
                                             ,
                                             Product([Symbol("Jinv_11", 1, GEO), Symbol("FE0_C1_D01_ip_j", 1, BASIS)])
                                            ])
                                        ,
                                        Symbol("", 2, CONST)
                                        ,
                                        Sum([
                                             Product([Symbol("Jinv_01", 1, GEO), Symbol("FE0_C1_D10_ip_k", 1, BASIS)])
                                             ,
                                             Product([Symbol("Jinv_11", 1, GEO), Symbol("FE0_C1_D01_ip_k", 1, BASIS)])
                                            ])
                                        ,
                                        Symbol("", 2, CONST)
                                       ])
                               ,
                               Product([
                                        Sum([
                                             Sum([
                                                  Product([Symbol("Jinv_01", 1, GEO), Symbol("FE0_C0_D10_ip_j", 1, BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_11", 1, GEO), Symbol("FE0_C0_D01_ip_j", 1, BASIS)])
                                                 ])
                                             ,
                                             Sum([
                                                  Product([Symbol("Jinv_00", 1, GEO), Symbol("FE0_C1_D10_ip_j", 1, BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_10", 1, GEO), Symbol("FE0_C1_D01_ip_j", 1, BASIS)])
                                                 ])
                                            ])
                                        ,
                                        Sum([
                                             Sum([
                                                  Product([Symbol("Jinv_01", 1, GEO), Symbol("FE0_C0_D10_ip_k", 1, BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_11", 1, GEO), Symbol("FE0_C0_D01_ip_k", 1, BASIS)])
                                                 ])
                                             ,
                                             Sum([
                                                  Product([Symbol("Jinv_00", 1, GEO), Symbol("FE0_C1_D10_ip_k", 1, BASIS)])
                                                  ,
                                                  Product([Symbol("Jinv_10", 1, GEO), Symbol("FE0_C1_D01_ip_k", 1, BASIS)])
                                                 ])
                                            ])
                                       ])
                              ])
                     ])
                     ,
                     Symbol("", 0.25, CONST)
                     ,
                     Symbol("W4_ip", 1, IP)
                     ,
                     Symbol("det", 1, GEO)
                     ])

        print
        start = time.time()
        expr_rem = expr.remove_nested()
        print "Elasticity2D: time, remove_nested():                  ", time.time() - start

        start = time.time()
        expr_exp = expr.expand()
        print "Elasticity2D: time, expand() (incl. remove_nested()): ", time.time() - start

        start = time.time()
        elasticity_exp = expand_operations(elasticity, get_format())
        print "Elasticity2D: time, old expand():                     ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "Elasticity2D: time, reduce_ops():                     ", time.time() - start

        start = time.time()
        elasticity_red = reduce_operations(elasticity, get_format())
        print "Elasticity2D: time, old reduce():                     ", time.time() - start

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

        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_rem)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(elasticity)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(elasticity_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(elasticity_red)))



    def testElasticityTerm(self):

        # expr:  0.25*W1*det*(FE0_C2_D001[0][j]*FE0_C2_D001[0][k]*Jinv_00*Jinv_21 + FE0_C2_D001[0][j]*FE0_C2_D001[0][k]*Jinv_00*Jinv_21)
        expr = Product([
                         Symbol('W1', 0.25, GEO), Symbol('det', 1, GEO),
                         Sum([Product([Symbol('FE0_C2_D001_0_j', 1, BASIS), Symbol('FE0_C2_D001_0_k', 1, BASIS),
                                       Symbol('Jinv_00', 1, GEO), Symbol('Jinv_21', 1, GEO)]),
                              Product([Symbol('FE0_C2_D001_0_j', 1, BASIS), Symbol('FE0_C2_D001_0_k', 1, BASIS),
                                  Symbol('Jinv_00', 1, GEO), Symbol('Jinv_21', 1, GEO)])
                             ])
                      ])

        print
        start = time.time()
        expr_rem = expr.remove_nested()
        print "ElasticityTerm: time, remove_nested():                  ", time.time() - start

        start = time.time()
        expr_exp = expr.expand()
        print "ElasticityTerm: time, expand() (incl. remove_nested()): ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "ElasticityTerm: time, reduce_ops():                     ", time.time() - start

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
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_rem)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))


    def testElasWeighted(self):
        expr = Product([
                          Symbol('W4', 1.0000, IP),
                          Symbol('det', 1.0000, GEO),
                          Sum([
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', 1.0000, BASIS),
                                        Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                        Symbol('Jinv_00', 1.0000, GEO),
                                        Symbol('Jinv_11', 1.0000, GEO),
                                        Symbol('w1', 1.0000, GEO)
                                        ]),
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', 1.0000, BASIS),
                                        Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                        Symbol('Jinv_01', 1.0000, GEO),
                                        Symbol('Jinv_10', 1.0000, GEO),
                                        Symbol('w0', 1.0000, GEO)
                                        ]),
                              Product([
                                        Symbol('w2', 1.0000, GEO),
                                        Sum([
                                              Product([
                                                      Symbol('FE0_C1_D01_ip_j', 1.0000, BASIS),
                                                      Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                                      Symbol('Jinv_00', 1.0000, GEO),
                                                      Symbol('Jinv_11', 1.0000, GEO),
                                                      Symbol('w1', 1.0000, GEO)
                                                      ]),
                                              Product([
                                                      Symbol('FE0_C1_D01_ip_j', 1.0000, BASIS),
                                                      Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                                      Symbol('Jinv_01', 1.0000, GEO),
                                                      Symbol('Jinv_10', 1.0000, GEO),
                                                      Symbol('w0', 1.0000, GEO)
                                                      ])
                                            ])
                                      ])
                              ])
                          ])
                                                       
        print
        start = time.time()
        expr_rem = expr.remove_nested()
        print "ElasWeighted: time, remove_nested():                  ", time.time() - start

        start = time.time()
        expr_exp = expr.expand()
        print "ElasWeighted: time, expand() (incl. remove_nested()): ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "ElasWeighted: time, reduce_ops():                     ", time.time() - start

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
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_rem)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))


    def testElasWeighted2(self):

        expr = Product([
                        Symbol('W4', 1.0000, IP),
                        Symbol('det', 1.0000, GEO),
                        Sum([
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', 1.0000, BASIS),
                                        Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                        Symbol('Jinv_00', 1.0000, GEO),
                                        Symbol('Jinv_10', 1.0000, GEO),
                                        Symbol('w1', 1.0000, GEO)
                                        ]),
                              Product([
                                        Symbol('FE0_C1_D01_ip_j', 1.0000, BASIS),
                                        Symbol('Jinv_01', 1.0000, GEO),
                                        Sum([
                                              Product([
                                                        Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                                        Symbol('Jinv_11', 1.0000, GEO),
                                                        Symbol('w0', 1.0000, GEO)
                                                        ]),
                                              Product([
                                                        Symbol('', 2.0000, CONST),
                                                        Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                                        Symbol('Jinv_11', 1.0000, GEO),
                                                        Symbol('w1', 1.0000, GEO)
                                                        ])
                                              ])
                                        ]),
                              Product([
                                        Symbol('w2', 1.0000, GEO),
                                        Sum([
                                            Product([
                                                    Symbol('FE0_C1_D01_ip_j', 1.0000, BASIS),
                                                    Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                                    Symbol('Jinv_00', 1.0000, GEO),
                                                    Symbol('Jinv_10', 1.0000, GEO),
                                                    Symbol('w1', 1.0000, GEO)
                                                    ]),
                                            Product([
                                                    Symbol('FE0_C1_D01_ip_j', 1.0000, BASIS),
                                                    Symbol('Jinv_01', 1.0000, GEO),
                                                    Sum([
                                                          Product([
                                                                  Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                                                  Symbol('Jinv_11', 1.0000, GEO),
                                                                  Symbol('w0', 1.0000, GEO)
                                                                  ]),
                                                          Product([
                                                                  Symbol('', 2.0000, CONST),
                                                                  Symbol('FE0_C1_D01_ip_k', 1.0000, BASIS),
                                                                  Symbol('Jinv_11', 1.0000, GEO),
                                                                  Symbol('w1', 1.0000, GEO)
                                                                  ])
                                                          ])
                                                    ])
                                            ])
                                        ])
                              ])
                        ])

        print
        start = time.time()
        expr_rem = expr.remove_nested()
        print "ElasWeighted2: time, remove_nested():                  ", time.time() - start

        start = time.time()
        expr_exp = expr.expand()
        print "ElasWeighted2: time, expand() (incl. remove_nested()): ", time.time() - start

        start = time.time()
        expr_red = expr_exp.reduce_ops()
        print "ElasWeighted2: time, reduce_ops():                     ", time.time() - start

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
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_rem)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_exp)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(expr_red)))
        self.assertAlmostEqual(eval(str(expr)), eval(str(opt_code)))


def suite():
    suite = unittest.TestSuite()
    # Classes and member functions
    suite.addTest(Tests('testSymbol'))
    suite.addTest(Tests('testProduct'))
    suite.addTest(Tests('testSum'))
    suite.addTest(Tests('testFraction'))
    suite.addTest(Tests('testMixedSymbols'))
    suite.addTest(Tests('testRemoveNested'))
    suite.addTest(Tests('testExpandOperations'))
    suite.addTest(Tests('testReduceVarType'))
    suite.addTest(Tests('testReduceOperations'))

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


#    print format
#    print BASIS
    if format == None:
#        print "none"
#        format = Format(FFC_OPTIONS).format
        set_format(Format(FFC_OPTIONS).format)

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(suite())

