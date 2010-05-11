"""This script runs a benchmark study on the form files found in the
current directory. It relies on the regression test script for
timings."""

__author__ = "Anders Logg"
__date__ = "2010-05-11"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

import os, glob

# Get list of test cases
test_cases = [f.split(".")[0] for f in glob.glob("*.ufl")]

# Run benchmark
os.chdir("../test/regression")
#os.system("python test.py --bench")

# Collect results
for test_case in test_cases:
    output = open("output/%s.out" % test_case).read()
    lines = [line for line in output.split("\n") if "bench" in line]
    if not len(lines) == 1:
        raise RuntimeError, "Unable to extract benchmark data from output"
    timing = float(lines[0].split(":")[-1])
    print test_case, timing
