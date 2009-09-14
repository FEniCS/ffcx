"Regression tests for FFC"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-05 -- 2009-08-21"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

import sys, shutil
from os import chdir, listdir, system, path, pardir, curdir, mkdir
from difflib import unified_diff

# Modified by Marie Rognes (meg), 2009-07-05
# Modified by Kristian B. Oelgaard, 2009-08-21

def check_forms(form_files, representation, exceptions):

    num_forms = len(form_files)
    forms_not_ok = []

    for form_file in form_files:
        if form_file in exceptions:
            continue

        print "Compiling and verifying form %s using %s representation ..." % (form_file, representation)
        # Compile form
        if system("python %s -fprecision=9 -s -r %s %s" % (path.join(pardir, "scripts", "ffc"), representation, form_file)) == 0:
            # Compare against reference
            code_file = form_file.split(".")[0] + ".h"
            f0 = open(path.join(pardir, "test", "regression", "reference", representation, code_file), "r")
            f1 = open(code_file, "r")
            c0 = f0.read().split("\n")
            c1 = f1.read().split("\n")
            f0.close()
            f1.close()
            if c0 != c1:
                print "*** Generated code does not match reference, diff follows"
                diff = unified_diff(c0, c1)
                for line in diff:
                    print line
                forms_not_ok += [form_file]

                # Copy failed file to temporary directory
                f_failed = open(path.join(pardir, "test", "regression", "tmp", representation, code_file), "w")
                f_failed.write("\n".join(c1))
                f_failed.close()
        else:
            forms_not_ok += [form_file]
    return forms_not_ok

# Check arguments, -nw specifices that we run in terminal mode
if len(sys.argv) == 2 and sys.argv[1] == "-nw":
    nw = True
else:
    nw = False

# Check both representations
representations = ["quadrature", "tensor"]

# Create temporary directory for generated files that fails
if not path.isdir("tmp"):
    mkdir("tmp")
    chdir("tmp")
    mkdir("quadrature")
    mkdir("tensor")
    chdir(pardir)

# Check all in demo directory
chdir(path.join(pardir, pardir, "demo"))
form_files = [f for f in listdir(curdir) if f.endswith(".ufl")]
form_files.sort()

# Some exceptions for the tensor representation
exceptions = {"tensor": ["Biharmonic.ufl",
                         "FunctionOperators.ufl",
                         "PoissonDG.ufl",
                         "QuadratureElement.ufl",
                         "TensorWeightedPoisson.ufl",
                         "Normals.ufl"],
              "quadrature" : []}

# Run regression tests for each representation
forms_not_ok = {}
for representation in representations:
    print "Checking %s representation \n" % representation
    forms_not_ok[representation] = check_forms(form_files, representation, exceptions[representation])
    print ""

# Print summary

test_failed = False
for representation in representations:
    num_forms = len(form_files) - len(exceptions[representation])
    num_forms_not_ok = len(forms_not_ok[representation])
    print "Summary for %s representation:" % representation
    if num_forms_not_ok == 0:
        print "\tRegression test completed for all %d forms ... All OK!" % num_forms
        print ""
    else:
        test_failed = True
        print "\t*** Regression test failed for %d of %d forms:" % (num_forms_not_ok, num_forms)
        for form_file in forms_not_ok[representation]:
            print "  " + form_file

        # Generate script for viewing diffs in meld                                    \
        file = open("../test/regression/viewdiff_%s.sh" % representation, "w")
        file.write("#!/bin/sh\n\n")
        for form_file in forms_not_ok[representation]:
            code_file = form_file.split(".")[0] + ".h"
            if nw:
                file.write("diff reference/%s/%s tmp/%s/%s\n" % (representation, code_file, representation, code_file))
            else:
                file.write("meld reference/%s/%s tmp/%s/%s\n" % (representation, code_file, representation, code_file))
        file.close()
        print "\nTo view diffs with meld, run the script viewdiff_%s.sh" % representation

#  Remove temporary directory if not errors occurred
if not test_failed:
    chdir(path.join(pardir, "test", "regression"))
    shutil.rmtree("tmp")

# Return error code if tests failed
sys.exit(test_failed)
