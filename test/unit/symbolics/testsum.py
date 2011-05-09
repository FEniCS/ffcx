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

class TestSum(unittest.TestCase):

    def testSum(self):
        "Test simple sum instance."

        f_0 = format["float"](0)
        f_1 = format["float"](1)
        f_2 = format["float"](2)
        f_3 = format["float"](3)

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

        self.assertEqual(repr(S0), "Sum([FloatValue(%s)])" % f_0)
        self.assertEqual(S0.t, CONST)
        self.assertEqual(repr(S1), "Sum([Symbol('x', BASIS)])")
#        self.assertEqual(repr(S4), "Sum([Symbol('x', BASIS), FloatValue(-2)])")
        self.assertEqual(repr(S4), "Sum([FloatValue(-%s), Symbol('x', BASIS)])" % f_2)
        self.assertEqual(repr(S9), "Sum([Symbol('y', GEO)])")

        self.assertEqual(str(S2), "(x + y)")
        self.assertEqual(str(S3), "(x + x)")
        self.assertEqual(str(S5), "(x + x-%s)" % f_2)
        self.assertEqual(str(S6), "(%s + x + x)" % f_1)
        self.assertEqual(str(S7), "(x + y-%s)" % f_2)
        self.assertEqual(str(S8), "(%s + x + x)" % f_3)
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

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestSum('testSum'))

