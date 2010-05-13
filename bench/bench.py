"""This script runs a benchmark study on the form files found in the
current directory. It relies on the regression test script for
timings."""

__author__ = "Anders Logg"
__date__ = "2010-05-11"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

import os, glob
from utils import print_table

# Test options
test_options = ["-r tensor", "-r tensor -O", "-r quadrature", "-r quadrature -O"]

# Get list of test cases
test_cases = sorted([f.split(".")[0] for f in glob.glob("*.ufl")])

# Open logfile
logfile = open("bench.log", "w")

# Iterate over options
os.chdir("../test/regression")
table = {}
for (j, test_option) in enumerate(test_options):

    # Run benchmark
    print "\nUsing options %s\n" % test_option
    os.system("python test.py --bench %s" % test_option)

    # Collect results
    for (i, test_case) in enumerate(test_cases):
        output = open("output/%s.out" % test_case).read()
        lines = [line for line in output.split("\n") if "bench" in line]
        if not len(lines) == 1:
            raise RuntimeError, "Unable to extract benchmark data for test case %s" % test_case
        timing = float(lines[0].split(":")[-1])
        table[(i, j)] = (test_case, test_option, timing)
        logfile.write("%s, %s, %g\n" % (test_case, test_option, timing))

# Close logfile
logfile.close()

# Print results
print_table(table, "FFC bench")
