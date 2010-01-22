#!/usr/bin/env python

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-01-06"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-06

# Pyhton modules
import unittest
import time

# FFC modules
from ffc.quadrature.symbolics import *
from ffc.quadrature.sumobj import _group_fractions
from ffc.cpp import format

class TestNotFinished(unittest.TestCase):

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
        e3 = _group_fractions(S3)
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

    if format == None:
        set_format(format)

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestNotFinished('testNotFinished'))

