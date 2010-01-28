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

