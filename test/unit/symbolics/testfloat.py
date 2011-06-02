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
set_float_formatting(FFC_PARAMETERS["precision"])

class TestFloat(unittest.TestCase):

    def testFloat(self):
            "Test simple FloatValue instance."
            f0 = FloatValue(1.5)
            f1 = FloatValue(-5)
            f2 = FloatValue(-1e-14)
            f3 = FloatValue(-1e-11)
            f4 = FloatValue(1.5)

    #        print "\nTesting FloatValue"
    #        print "f0: '%s'" %f0
    #        print "f1: '%s'" %f1
    #        print "f2: '%s'" %f2
    #        print "f3: '%s'" %f3

            self.assertEqual(repr(f0), "FloatValue(%s)" % format["float"](1.5))
            self.assertEqual(repr(f1), "FloatValue(%s)" % format["float"](-5))
            self.assertEqual(repr(f2), "FloatValue(%s)" % format["float"](0))
            self.assertEqual(repr(f3), "FloatValue(%s)" % format["float"](-1e-11))

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

if __name__ == "__main__":

    # Run all returned tests
    runner = unittest.TextTestRunner()
    runner.run(TestFloat('testFloat'))
