"Regression tests for FFC"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-05 -- 2007-03-05"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

from os import chdir, listdir, system
from difflib import unified_diff

# Check all in demo directory
chdir("../../demo")
form_files = listdir(".")
# Temporary until all forms work
form_files = ["Poisson.form", "Mass.form"]
num_forms = len(form_files)
num_forms_ok = 0

# Iterate over form files
for form_file in form_files:

    # Check only files ending in .form
    if ".form" in form_file:
        print "Compiling and verifying form %s..." % form_file

        # Compile form
        if system("../bin/ffc -s -l ufc %s" % form_file) == 0:

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

# Print summary
if num_forms_ok == num_forms:
    print ""
    print "Regression test completed for all %d forms" % num_forms_ok
    print ""
    print "OK"
else:
    print ""
    print "Regression test failed for %d of %d forms" % (num_forms - num_forms_ok, num_forms)
