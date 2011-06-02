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

class TestSymbol(unittest.TestCase):

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

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestSymbol('testSymbol'))

