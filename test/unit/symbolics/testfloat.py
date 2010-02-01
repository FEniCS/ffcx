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
from ffc.cpp import format, set_float_formatting
from ffc import default_parameters
set_float_formatting(default_parameters()['precision'])

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
