#!/usr/bin/env python

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-01-06"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-02-01

# Pyhton modules
import unittest
import time

# FFC modules
from ffc.quadrature.symbolics import *
from ffc.quadrature.sumobj import _group_fractions
from ffc.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])

class TestNotFinished(unittest.TestCase):

    def testNotFinished(self):
        "Stuff that would be nice to implement."

        f_1 = format["float"](1)
        f_2 = format["float"](2)
        f_4 = format["float"](4)
        f_8 = format["float"](8)

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
        e3 = _group_fractions(S3)
        e4 = Sum([Fraction(f1*s0, a*b*c), Fraction(s0, a*b)]).expand().reduce_ops()

        # Tests that pass the current implementation
        self.assertEqual(str(e0), '%s/(%s*x + %s*y)' % (f_4, f_2, f_8))
        self.assertEqual(str(e1), 'x/(x + x*y)')
        self.assertEqual(str(e2), '(%s + y)/(x + x*y)' % f_1)
        self.assertEqual(str(e3), '%s/x' % f_8)
        self.assertEqual(str(e4), 'x*(%s/(a*b) + %s/(a*b*c))' % (f_1, f_2))

        # Tests that should pass in future implementations (change NotEqual to Equal)
        self.assertNotEqual(str(e0), '%s/(x + %s*y)' % (f_2, f_4))
        self.assertNotEqual(str(e1), '%s/(%s + y)' % (f_1, f_1))
        self.assertNotEqual(str(e2), '%s/x' % f_1)
        self.assertNotEqual(str(e4), 'x*(%s/c + %s)/(a*b)' % (f_2, f_1))

        # TODO: Would it be a good idea to reduce expressions wrt. var_type
        # without first expanding?
        E0 = Product([ Sum([ Product([ Symbol('B0', BASIS), Product([Symbol('B1', BASIS), Sum([s0]), Sum([s0])]) ]),
                             Product([Symbol('B0', BASIS), Symbol('B1', BASIS)]) ]) ])
        Er0 = E0.reduce_vartype(BASIS)
        Ex0 = E0.expand().reduce_vartype(BASIS)
#        print "%s, red(BASIS): ('%s', '%s')" %(E0, Ex0[0][0], Ex0[0][1])
#        print "%s, red(BASIS): ('%s', '%s')" %(E0, Er0[0], Er0[1])
        self.assertNotEqual( Ex0[0][1], Er0[1].expand() )

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestNotFinished('testNotFinished'))

