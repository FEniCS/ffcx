#!/usr/bin/env python

"Test suite for FFC: Verify tensors tabulated by forms in the demo directory"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2008-07-11 -- 2008-09-10"
__copyright__ = "Copyright (C) 2008 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

import sys
import getopt
from glob import glob
from os import chdir, getcwd, path, system
import ufc_benchmark
import numpy
import pickle
import math
import time

from ffc.jit.jit import jit
# Enable ffc without installation
#sys.path.append("../../")
#from ffc import *

def main(argv):
    "Main function, handle arguments and run tests accordingly."

    # Get command-line arguments
    try:
        opts, args = getopt.getopt(argv, "hnr:t:T:", \
        ["help", "new_references", "representation=", "tolerance=", "type="])
    except getopt.GetoptError:
        usage()
        return 2

    # Run tests for both representations as default
    representations = ["quadrature", "tensor"]

    # Default options
    tolerance = 1e-14
    new_references = False
    form_types = ["form", "ufl", "all"]
    form_type = "form"

    # Get options
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            return 0
        elif opt in  ("-r", "--representation"):
            if arg in representations:
                representations = [arg]
            else:
                usage()
                return 2
        elif opt in  ("-n", "--new_references"):
            new_references = True
            representations = ["tensor"]
        elif opt in  ("-t", "--tolerance"):
            try:
                tolerance = float(arg)
            except: usage()
        elif opt in  ("-T", "--type"):
            if arg in form_types:
                form_type = arg
            else:
                usage()
                return 2
        else:
            usage()
            return 2

    test_options = {"tolerance": tolerance,
                    "new_references": new_references,
                    "form_files": args,
                    "form_type": form_type}

    # Print test options
    print "\nThe following test options will be used"
    print "====================================================================\n"
    for key, val in test_options.items():
        if key == "form_files" and not val:
            print (key, "all")
        else:
            print (key, val)
    print "representations: ", representations
    print "\n====================================================================\n"

    # Clean test_cache
    print "\nDeleting forms in cache"
    print "====================================================================\n"
    system("rm -rf form_*")
    system("rm -rf test_cache/form_*")
    print "Done"
    print "\n====================================================================\n"

    # Run tests and get summary
    summary = {}
    for representation in representations:
        test_options["representation"] = representation
        summary[representation] = run_tests(test_options)

    # Print summary for each test
    for representation in representations:
        print "\n*** Summary for %s representation ***" % representation
        print "====================================================================\n"
        print summary[representation]
        print "====================================================================\n"
    print ""

    return 0

def run_tests(test_options):
    "Run tests for given options."

    new_references = test_options["new_references"]

    if not new_references:
        print "\nVerifying tensors using %s representation" % test_options["representation"]
        print "Tolerance for norm when comparing tensors: %.g" % test_options["tolerance"]
        print "====================================================================\n"
    else:
        print "\n*** Generating new reference tensors using %s representation ***" % test_options["representation"]
        print "====================================================================\n"

    # Get form files from user and demo directory
    form_files = test_options["form_files"]
    chdir("../../demo")
    form_type = test_options["form_type"]
    demo_files = []
    if form_type == "form" or form_type == "all":
        demo_files += glob("*.form")
    if form_type == "ufl" or form_type == "all":
        demo_files += glob("*.ufl")
    chdir("../test/verify_tensors")

    # If not form files are specified, check all in demo directory
    if not form_files:
        print "Verifying tensors for all forms in the demo directory"
        form_files = demo_files

    # Check that all form files are present in demo directory and remove forms
    # that are known to break the test
    do_not_compile = ["TensorWeightedPoisson.ufl", "MixedPoisson.ufl", "VectorLaplaceGradCurl.ufl"]

    working_forms = ["Constant.ufl", "Elasticity.ufl", "EnergyNorm.ufl", "Equation.ufl", "FunctionOperators.ufl",
                     "Heat.ufl", "Mass.ufl", "NavierStokes.ufl", "NeumannProblem.ufl", "Optimization.ufl",
                     "Poisson.ufl", "PoissonSystem.ufl", "SubDomain.ufl", "SubDomains.ufl",
                     "PoissonDG.ufl"]

    all_forms = ["Constant.ufl", "Mass.ufl", "P5tet.ufl", "QuadratureElement.ufl",
                 "Elasticity.ufl", "MixedMixedElement.ufl", "P5tri.ufl", "Stokes.ufl",
                 "EnergyNorm.ufl", "MixedPoisson.ufl", "PoissonDG.ufl", "SubDomains.ufl",
                 "Equation.ufl", "NavierStokes.ufl", "PoissonSystem.ufl", "SubDomain.ufl",
                 "FunctionOperators.ufl", "NeumannProblem.ufl", "Poisson.ufl", "TensorWeightedPoisson.ufl",
                 "Heat.ufl", "Optimization.ufl", "Projection.ufl", "VectorLaplaceGradCurl.ufl"]

    new_files = []
    for form_file in form_files:
        if form_file in demo_files and form_file not in new_files:
            if form_file not in do_not_compile:
                new_files.append(form_file)
            else:
                print "*** %s is known to break the test, omitting from test." % form_file
        else:
            do_not_compile.append(form_file)
            print "*** %s is not in demo directory, omitting from test." % form_file

    form_files = new_files
    form_files.sort()

    # Get total number of forms
    num_forms = len(form_files)

    # Number of forms that involve tabulating tensors (will be reduced accordingly)
    num_tensor_forms = num_forms

    # Forms that pass everything
    num_forms_ok = 0

    # Forms that are not processed correctly (FFC syntax)
    forms_not_read_ok = []

    # Forms that are not compiled correctly by FFC (jit)
    forms_not_compiled_ok = []

    # Forms for which the norm of the tensor compared to reference was too large
    forms_not_compared_ok = []

    # Iterate over form files
    for form_file in form_files:

        # Get suffix and set compiler type
        suffix = form_file.split(".")[-1]
        if suffix == "ufl":
            test_options["compiler"] = "ufl"
        elif suffix == "form":
            test_options["compiler"] = "ffc"

        print ""
        if not new_references:
            print "\nCompiling and verifying tensors for %s..." % form_file
            print "--------------------------------------------------------------------"
        else:
            print "\nGenerating new tensors for %s..." % form_file
            print "--------------------------------------------------------------------"

        # Read the forms from the form file
        (forms, read_ok) = read_forms("../../demo", form_file)
        if not read_ok:
            forms_not_read_ok.append(form_file)
        else:
            # No forms were present, reduce number of tensor forms
            if not forms:
                print "\nThis form does not specify any integral, omitting from test."
                num_tensor_forms -= 1
            else:
                # Verify correctness of forms
                (ok_compile, norm_ok) = verify_forms(form_file, forms, forms_not_compiled_ok, forms_not_compared_ok, test_options)
                # Increase number of correct forms if all tests were OK
                if norm_ok and ok_compile:
                    num_forms_ok += 1

    print "\n\nFinished tests for %s representation" % test_options["representation"]
    print "\n====================================================================\n"

    # Generate summary
    summary = "\n"

    # Some forms might not involve tabulating tensors (only generate element)
    if num_forms != num_tensor_forms:
        summary += "%d out of %d forms involved tabulating tensors\n" % (num_tensor_forms, num_forms)

    # Print list of forms that we're omitted from the tests
    if do_not_compile:
        summary += "The following forms were not included in the test either because they are"
        summary += "\nknown to break the test or they are not present in the ../../demo directory:\n"
        summary += "\n".join(do_not_compile)
        summary += "\n\n"

    # If all tests succeded there's no need to add detailed info
    if num_forms_ok == num_tensor_forms:
        # If we ran tests
        if not test_options["new_references"]:
            summary += "Test completed for all %d forms\n\nOK\n" % num_forms_ok
        else:
            summary += "Generated new references for all %d forms\n\nOK\n" % num_forms_ok
    else:
        # If there was problems reading some forms (interpreting FFC syntax)
        # this is highly unlikely as it would be picked up by the other tests
        # but it's a good test if this module is used elsewhere
        if forms_not_read_ok:
            summary += "The following forms do not have correct FFC syntax:\n"
            summary += "\n".join(forms_not_read_ok)
            summary += "\n\n"

        # If some forms could not be compiled by FFC report them. This is
        # particullary relevant for quadrature representation
        if forms_not_compiled_ok:
            summary += "The following forms could not be compiled by FFC:\n"
            summary += "\n".join(forms_not_compiled_ok)
            summary += "\n\n"

        # If some forms resulted in too large differences compared to the reference
        # solutions, report.
        if forms_not_compared_ok:
            max_len = max([len(f) for f,v in forms_not_compared_ok])
            summary += "The norm of the tensor(s) compared to reference was too large for the following forms:\n"
            for form, norm in forms_not_compared_ok:
                space = max_len - len(form)
                summary += form + " "*space + " norm: %.4e, tolerance: %.4e" % (norm, test_options["tolerance"]) + "\n"
            summary += "\n"

    return summary

def read_forms(demo_dir, form_file):
    "Read/construct forms in form file"

    # Create the correct type of script
    # Get filename suffix
    suffix = form_file.split(".")[-1]
    form_file = form_file.split(".")[0]
    script = ""
    # Check file suffix and parse file/generate module
    if suffix == "ufl":
        script = _make_script_ufl(path.join(demo_dir, form_file))
    elif suffix == "form":
        script = _make_script_form(path.join(demo_dir, form_file))

    # Initialise bilinear and linear forms
    a, L, M = (0, 0, 0)
    forms = {}
    read_ok = 0
    try:
        # FIXME: there must be a cleverer way of doing this....
        global_dict = globals()
        global_dict['a']= a
        global_dict['L']= L
        global_dict['M']= M
        print ""
        execfile(script, global_dict)
        for f in ['a','L','M']:
            if global_dict[f] != 0:
                forms[f] = global_dict[f]
        read_ok = 1
    except Exception, what:
        print "*** An error occured while reading form file"
        # Reset forms
        forms = {}
        print "What: ", what
        pass

    return (forms, read_ok)

# New version for .ufl files
def _make_script_ufl(filename):
    "Create Python script from given .ufl file and return name of script"

    print "Generating Python script for ", filename + ".ufl"
    # Read input
    infile = open(filename + ".ufl", "r")
    input = infile.read()
    infile.close()

    script = filename + ".py"

    # Generate output
    output = """from ufl import *
%s
""" % input

    # Write output
    outfile = open(script, "w")
    outfile.write(output)
    outfile.close()

    # Return script filename
    return script

# New version for .ufl files
def _make_script_form(filename):
    "Create Python script from given .ufl file and return name of script"

    print "Generating Python script for ", filename + ".form"
    # Read input
    infile = open(filename + ".form", "r")
    input = infile.read()
    infile.close()

    script = filename + ".py"

    # Generate output
    output = """from ffc import *
%s
""" % input

    # Write output
    outfile = open(script, "w")
    outfile.write(output)
    outfile.close()

    # Return script filename
    return script

def verify_forms(form_file, forms, forms_not_compiled_ok, forms_not_compared_ok, test_options):
    # Some flags that indicates success/failure
    ok_compile = True
    norm_ok = True
    for form_type, form in forms.items():
        (ok_compile, norm_ok) = verify_form(form_file, form_type, form,\
          forms_not_compiled_ok, forms_not_compared_ok, ok_compile, norm_ok, test_options)
    return (ok_compile, norm_ok)

def verify_form(form_file, form_type, form, forms_not_compiled_ok, forms_not_compared_ok, ok_compile, norm_ok, test_options):

    compiled_form, module, form_data = (0, 0, 0)

    try:
        # Compile the form with jit
        opt = {"representation":test_options["representation"], "cache_dir":"test_cache", "compiler": test_options["compiler"]}
        (compiled_form, module, form_data) = jit(form, opt)
    except Exception, what:
        ok_compile = False
        forms_not_compiled_ok.append(form_file)
        print "\n*** An error occured while compiling form"
        print "What: ", what
        pass

    norm = 1.0
    try:
        # Compute norm of tensor compared to reference and compare to tolerance
        ref_file = path.splitext(form_file)[0] + "_" + form_type
        norm = compute_norm(compiled_form, form_data, ref_file, test_options)
        print "\nNorm for %s, %s: " %(form_file, form_type) , norm
        if norm > test_options["tolerance"]:
            forms_not_compared_ok.append((form_file + ", " + form_type + ":", norm))
            norm_ok = False
    except Exception, what:
        norm_ok = False
        forms_not_compared_ok.append((form_file, norm))
        print "An error occured while computing norm"
        print "What: ", what
        pass

    return (ok_compile, norm_ok)

def compute_norm(compiled_form, form_data, file_name, test_options):

    # Initialise variables
    norm = 0.0
    num_coefficients, cell_shape, mesh, cell = (0, 0, 0, 0)

    rank = form_data.rank
    if test_options["compiler"] == "ffc":
        cell_shape = form_data.cell_dimension
        num_coefficients = form_data.num_coefficients
    else:
        cell_shape = form_data.geometric_dimension
        num_coefficients = form_data.num_functions
    num_arguments = num_coefficients + rank

    # FIXME: Can only handle, interval, triangle, tetrahedron
    # Random cells were generated by:
    # >>> random.uniform(-2.5,2.5)
    if cell_shape == 1:
        mesh = ufc_benchmark.Mesh(1, 1, [2, 1, 0])
        # Reference cell
        # cell = ufc_benchmark.Cell(1, 1, [[0], [1]], [2, 1, 0, 0])
        # Random cell
        cell = ufc_benchmark.Cell(1, 1, [[-1.445], [0.4713]], [2, 1, 0, 0])
    elif cell_shape == 2:
        mesh = ufc_benchmark.Mesh(2, 2, [3, 3, 1])
        # Reference cell
        # cell = ufc_benchmark.Cell(2, 2, [[0, 0], [1, 0], [1, 1]], [3, 3, 1, 0])
        # Random cell
        cell = ufc_benchmark.Cell(2, 2, [[-2.2304, -0.88317], [1.3138, -1.0164],\
                                         [0.24622, 1.4431]], [3, 3, 1, 0])
    elif cell_shape == 3:
        mesh = ufc_benchmark.Mesh(3, 3, [4, 6, 4, 1])
        # Reference cell
        # cell = ufc_benchmark.Cell(3, 3, [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [4, 6, 4, 1])
        # Random cell
        cell = ufc_benchmark.Cell(3, 3, [[-2.2561, -1.6144, -1.7349], [-1.5612, -1.5121, -0.17675],\
                                         [1.6861, -1.1494, 2.4070], [0.52083, 1.1940, 1.8220]], [4, 6, 4, 1])
    else:
        print "*** Unknown cell dimension %d, returning norm = 1.0" % cell_shape
        return 1.0

    # Initialise dofmaps
    dof_maps = [0]*num_arguments
    for i in range(num_arguments):
        dof_maps[i] = compiled_form.create_dof_map(i)
        dof_maps[i].init_mesh(mesh)

    # Get number of integrals for given form
    num_cell_integrals = compiled_form.num_cell_integrals()
    num_exterior_facet_integrals = compiled_form.num_exterior_facet_integrals()
    num_interior_facet_integrals = compiled_form.num_interior_facet_integrals()

    # Simple coefficient generator
    w = [0]*num_coefficients
    for i in range(num_coefficients):
        w[i] = [0]*(dof_maps[rank+i].local_dimension())
        for j in range(dof_maps[rank+i].local_dimension()):
            w[i][j] = 1.111 + (i + j)/1.111
    macro_w = [0]*num_coefficients
    for i in range(num_coefficients):
        macro_w[i] = [0]*(2*dof_maps[rank+i].local_dimension())
        for j in range(2*dof_maps[rank+i].local_dimension()):
            macro_w[i][j] = 1.111 + (i + j)/1.111

    A = numpy.array([0.0])
    # Add contributions from ALL domains from cell integrals
    if num_cell_integrals:
        domain = 0
        # Get shape of A and reset values
        try:
            A = numpy.array(ufc_benchmark.tabulate_cell_integral(compiled_form, w, cell, domain))
            A = numpy.zeros(numpy.shape(A))
            for domain in range(num_cell_integrals):
                A += ufc_benchmark.tabulate_cell_integral(compiled_form, w, cell, domain)
        except Exception, what:
            print "*** An error occured while calling tabulate_foo_integral(), returning norm = 1.0"
            print "What: ", what
            return 1.0

    # Add contributions from ALL domains and facets from exterior integrals
    if num_exterior_facet_integrals:
        domain, facet = (0, 0)
        try:
            if not A.any():
                A = numpy.array(ufc_benchmark.tabulate_exterior_facet_integral(compiled_form, w, cell, facet, domain))
                A = numpy.zeros(numpy.shape(A))
            for domain in range(num_exterior_facet_integrals):
                for facet in range(cell.num_entities[cell_shape - 1]):
                    A += ufc_benchmark.tabulate_exterior_facet_integral(compiled_form, w, cell, facet, domain)
        except Exception, what:
            print "*** An error occured while calling tabulate_foo_integral(), returning norm = 1.0"
            print "What: ", what
            return 1.0

    # Add contributions from ALL domains and facets from interior integrals
    # FIXME: this currently makes no sense (integrating interior facets on 1 cell)
    #        but it should be OK since we just compare numbers.
    macro_A = numpy.array([0.0])
    if num_interior_facet_integrals:
        domain, facet0, facet1 = (0,0,0)
        try:
            macro_A = numpy.array(ufc_benchmark.tabulate_interior_facet_integral(compiled_form, macro_w, cell, cell, facet0, facet1, domain))
            macro_A = numpy.zeros(numpy.shape(macro_A))
            for domain in range(num_interior_facet_integrals):
                for facet in range(cell.num_entities[cell_shape - 1]):
                    macro_A += ufc_benchmark.tabulate_interior_facet_integral(compiled_form, macro_w, cell, cell, facet0, facet1, domain)
        except Exception, what:
            print "*** An error occured while calling tabulate_foo_integral(), returning norm = 1.0"
            print "What: ", what
            return 1.0

    # Add A to the upper left quadrant of macro_A, it makes no sense,
    # but the numbers should be OK
    if not macro_A.any():
        macro_A = A
    elif A.any():
        dims = numpy.shape(A)
        for i in range(dims[0]):
            for j in range(dims[1]):
                macro_A[i][j] += A[i][j]

    # We need to compare tensor to reference tensor
    if not test_options["new_references"]:
        A_ref = 0.0
        try:
            f = open(path.join("references", file_name), "r")
            A_ref = pickle.load(f)
            f.close()
        except Exception, what:
            print "*** An error occured while trying to load reference value: %s" % path.join("references", file_name)
            print "*** Maybe you need to generate the reference? Returning norm = 1.0"
            print "What: ", what
            return 1.0

        # Compute norm
        # Compute difference
#        print "macro_A: ", macro_A
#        print "A_ref: ", A_ref
        C = macro_A - A_ref

        # Normalise differences with reference values, sum tensor and take
        # absolute value
        s = numpy.sum(C/(1 + A_ref))
        norm = math.sqrt(s*s)

    # We need to generate a new reference tensor
    else:
        f = open(path.join("references", file_name), "w")
        if not macro_A.any():
            print "*** Warning: reference tensor only contains zeros!"
        pickle.dump(macro_A, f)
        f.close()

    return norm

def usage():
    "Display usage info."
    print """\nUsage: ./test.py [OPTION] [FILES (optional)], or: python test.py [OPTION] [FILES (optional)]

  -h, --help              Display this info.


  -n, --new_references    Generate new reference tensors for all forms in
                          [FILES] using tensor representation.

                          Example:

                            ./test.py -n Poisson.form

                          This generates new references for ../../demo/Poisson.form
                          and put the result in ./references

                          IMPORTANT! This option should only be used in one of
                          the following cases:

                          0. If a new form file was added to the ../../demo
                             directory.

                          1. If a bug was discovered, and the output to
                             tabulate_tensor() has changed.

                          2. If the benchmark has changed for some reason,
                             like the element which is being integrated.

                          3. If a form file in the ../../demo directory has changed.

                          4. If there's another good reason, put it here...


  -r, --representation    Specify the representation ('tensor' or 'quadrature') for
                          which to run the tests.

                          Example:

                            ./test.py -r quadrature

                          This will test ALL forms in ../../demo using quadrature representation.

                            ./test.py -r tensor Poisson.form

                          This will test ../../demo/Poisson.form using tensor representation.

                          If no representation is specified the tests will be
                          run for both representations.


  -t, --tolerance         Specify the tolerance that should be used when comparing
                          tensors.

                          Example:

                            ./test.py -t 1e-10

                          This will use a tolerance of 1e-10 when comparing the
                          norm between tensors being computed and those in the
                          ./references directory.

                          If no tolerance is specified, the default value will be
                          used (1e-14).


  -T, --type              FIXME: Remove this option once the .form format is removed!

                          Specify which type of forms to test. Can be either
                          the 'native' FFC i.e. 'form' (default), UFL i.e.
                          'ufl' or both types i.e. 'all'. E.g.,

                          ./test.py -T ufl
"""
    return

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
