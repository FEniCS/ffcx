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

class TestProduct(unittest.TestCase):

    def testProduct(self):
        "Test simple product instance."

        f_0 = format["float"](0)
        f_1 = format["float"](1)
        f0 = FloatValue(-2.0)
        f1 = FloatValue(3.0)
        f2 = FloatValue(0)
        f3 = FloatValue(-1)
        f4 = FloatValue(1)
        f5 = FloatValue(-0.5)
        f6 = FloatValue(2.0)
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
        p7 = Product([s0, f4, s1]).expand().reduce_ops()
        p8 = Product([s0, f0, s2, f5])
        p9 = Product([s0, s1])
        p10 = Product([p0, p1])
        p11 = Product([f5, f0])
        p12 = Product([f6, f5])
        p13 = Product([f6, f5]).expand()
        p14 = Product([f1, f2])
        p_tmp = Product([f1])
        p_tmp.expand()
        p15 = Product([p_tmp, s0])

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
#        print "\np9: %s * %s = '%s'" % (s0, s1, p9)
#        print "\np10: %s * %s = '%s'" % (p0, p1, p10)
#        print "\np11: %s * %s = '%s'" % (f6, f1, p11)
#        print "\np12: %s * %s = '%s'" % (f6, f5, p12)
#        print "\np13: %s * %s = '%s'" % (f6, f5, p13)
#        print "\np14: %s * %s = '%s'" % (f1, f2, p14)

        self.assertEqual(repr(p0), "Product([FloatValue(%s)])" % f_0)
        self.assertEqual(repr(p1), "Product([Symbol('x', BASIS)])")
        self.assertEqual(repr(p3), "Product([FloatValue(%s), Symbol('x', BASIS), Symbol('y', GEO)])"\
                                    % format["float"](3))
        self.assertEqual(repr(p6), "Product([FloatValue(-%s), Symbol('x', BASIS), Symbol('y', GEO)])" % f_1)
        self.assertEqual(repr(p7), "Product([Symbol('x', BASIS), Symbol('y', GEO)])")
        self.assertEqual(repr(p8), "Product([Symbol('x', BASIS), Symbol('z', GEO)])")
        self.assertEqual(str(p2), 'x*y')
        self.assertEqual(str(p4), '%s' % f_0)
        self.assertEqual(str(p5), '-%s*x*y*z' % format["float"](6))
        self.assertEqual(str(p6), ' - x*y')
        self.assertEqual(str(p7), 'x*y')
        self.assertEqual(str(p8), 'x*z')
        self.assertEqual(str(p9), 'x*y')
        self.assertEqual(p0.val, 0)
        self.assertEqual(str(p10), '%s' % f_0)
        self.assertEqual(str(p11), '%s' % f_1)
        self.assertEqual(str(p12), '-%s' % f_1)
        self.assertEqual(str(p13), '-%s' % f_1)
        self.assertEqual(repr(p14), "Product([FloatValue(%s)])" % f_0)
        self.assertEqual(repr(p14.expand()), "FloatValue(%s)" % f_0)

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
        self.assertEqual(p9.ops(), 1)
        self.assertEqual(p10.ops(), 0)
        self.assertEqual(p14.ops(), 0)

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

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestProduct('testProduct'))

