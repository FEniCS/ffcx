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
from ffc.log import push_level, pop_level, CRITICAL

class TestFraction(unittest.TestCase):

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

        # Silence output
        push_level(CRITICAL)
        self.assertRaises(Exception, Fraction, f0, f2)
        self.assertRaises(Exception, Fraction, s0, f2)
        pop_level()

        self.assertEqual(repr(F0), "Fraction(FloatValue(%s), FloatValue(%s))"\
                                    % (format["float"](-1.5), format["float"](1)))
        self.assertEqual(repr(F2), "Fraction(Symbol('x', BASIS), Symbol('y', GEO))")

        self.assertEqual(str(F0), "%s" % format["float"](-1.5))
        self.assertEqual(str(F1), "%s" % format["float"](0))
        self.assertEqual(str(F2), "x/y")
        self.assertEqual(str(F3), "x/%s" % format["float"](3))
        self.assertEqual(str(F4), "-%s/y" % format["float"](2))
        self.assertEqual(str(F5), "%s" % format["float"](0))

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

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestFraction('testFraction'))

