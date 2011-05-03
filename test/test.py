"""Run all tests, including unit tests and regression tests"""

# Copyright (C) 2007 Anders Logg
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
# First added:  2007-06-09
# Last changed: 2007-06-09

import os
import re
import sys

pwd = os.path.dirname(os.path.abspath(__file__))

# Tests to run
tests = ["unit", "regression"]

failed = []

# Run tests
for test in tests:
    print "Running tests: %s" % test
    print "----------------------------------------------------------------------"
    os.chdir(os.path.join(pwd, test))
    failure = os.system("python test.py")
    if failure:
        failed.append(test)
    print ""

sys.exit(len(failed))
