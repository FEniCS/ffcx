"Regression tests for FFC"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-05 -- 2008-04-29"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

import sys
from os import chdir, listdir, system
from difflib import unified_diff

# Check arguments, -nw specifices that we run in terminal mode
if len(sys.argv) == 2 and sys.argv[1] == "-nw":
    nw = True
else:
    nw = False

# Check all in demo directory
chdir("../../demo")
form_files = [f for f in listdir(".") if f[-5:] == ".form"]
form_files.sort()
num_forms = len(form_files)
num_forms_ok = 0
forms_not_ok = []

# Iterate over form files
for form_file in form_files:

    print "Compiling and verifying form %s..." % form_file

    # Compile form
    if system("../scripts/ffc -fprecision=9 -s %s" % form_file) == 0:

        # Compare against reference
        code_file = form_file.split(".")[0] + ".h"
        f0 = open("../test/regression/reference/%s" % code_file, "r")
        f1 = open(code_file, "r")
        c0 = f0.read().split("\n")
        c1 = f1.read().split("\n")
        f0.close()
        f1.close()
        if c0 == c1:
            num_forms_ok += 1
        else:
            print "*** Generated code does not match reference, diff follows"
            diff = unified_diff(c0, c1)
            for line in diff:
                print line
            forms_not_ok += [form_file]

# Print summary
if num_forms_ok == num_forms:
    print ""
    print "Regression test completed for all %d forms" % num_forms_ok
    print ""
    print "OK"
else:
    print ""
    print "*** Regression test failed for %d of %d forms:" % (num_forms - num_forms_ok, num_forms)
    for form_file in forms_not_ok:
        print "  " + form_file

    # Generate script for viewing diffs in meld                                    \
    file = open("../test/regression/viewdiff.sh", "w")
    file.write("#!/bin/sh\n\n")
    for form_file in forms_not_ok:
        code_file = form_file.split(".")[0] + ".h"
        if nw:
            file.write("diff reference/%s ../../demo/%s\n" % (code_file, code_file))
        else:
            file.write("meld reference/%s ../../demo/%s\n" % (code_file, code_file))
    file.close()
    print "\nTo view diffs with meld, run the script viewdiff.sh"

    # Return error code if tests failed
    sys.exit(True)
