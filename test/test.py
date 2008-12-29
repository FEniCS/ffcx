"""Run all tests, including unit tests and regression tests"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-06-09 -- 2007-06-09"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

import os
import re

pwd = os.path.dirname(os.path.abspath(__file__))

# Tests to run
tests = ["unit", "regression"]

# Run tests
for test in tests:
    print "Running tests: %s" % test
    print "----------------------------------------------------------------------"
    os.chdir(os.path.join(pwd, test))
    os.system("python test.py")
    print ""
