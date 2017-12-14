# -*- coding: utf-8 -*-
"""This script runs a benchmark study on the form files found in the
current directory. It relies on the regression test script for
timings."""

# Copyright (C) 2010 Anders Logg
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

import os
import glob
import sys
from utils import print_table

# Test options
test_options = ["-r tensor",
                #"-r tensor -O",
                "-r quadrature",
                "-r quadrature -O",
                "-r uflacs"]

# Get list of test cases
test_cases = sorted([f.split(".")[0] for f in glob.glob("*.ufl")])

# Open logfile
logfile = open("bench.log", "w")

# Iterate over options
os.chdir("../test/regression")
table = {}
for (j, test_option) in enumerate(test_options):

    # Run benchmark
    print("\nUsing options %s\n" % test_option)
    os.system(sys.executable + " test.py --bench %s" % test_option)

    # Collect results
    for (i, test_case) in enumerate(test_cases):
        output = open("output/%s.out" % test_case).read()
        lines = [line for line in output.split("\n") if "bench" in line]
        if not len(lines) == 1:
            raise RuntimeError("Unable to extract benchmark data for test case %s" % test_case)
        timing = float(lines[0].split(":")[-1])
        table[(i, j)] = (test_case, test_option, timing)
        logfile.write("%s, %s, %g\n" % (test_case, test_option, timing))

# Close logfile
logfile.close()

# Print results
print_table(table, "FFC bench")
