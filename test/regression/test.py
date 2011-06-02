"""This script compiles and verifies the output for all form files
found in the 'demo' directory. The verification is performed in two
steps. First, the generated code is compared with stored references.
Then, the output from all functions in the generated code is compared
with stored reference values.

This script can also be used for benchmarking tabulate_tensor for all
form files found in the 'bench' directory. To run benchmarks, use the
option -b.
"""

# Copyright (C) 2010 Anders Logg, Kristian B. Oelgaard and Marie E. Rognes
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
#
# First added:  2010-01-21
# Last changed: 2010-05-12

# FIXME: Need to add many more test cases. Quite a few DOLFIN
# FIXME: forms failed after the FFC tests passed. Also need to
# FIXME: to check with different compiler options, representation
# FIXME: and DOLFIN wrappers.

import os, sys, shutil, commands, difflib
from numpy import array, shape, abs, max
from ffc.log import begin, end, info, info_red, info_green, info_blue
from ufctest import generate_test_code

# Parameters
tolerance = 1e-9
output_directory = "output"
demo_directory = "../../../demo"
bench_directory = "../../../bench"

# Global log file
logfile = None

def run_command(command):
    "Run command and collect errors in log file."
    (status, output) = commands.getstatusoutput(command)
    if status == 0:
        return True
    global logfile
    if logfile is None:
        logfile = open("../error.log", "w")
    logfile.write(output + "\n")
    return False

def log_error(message):
    "Log error message."
    global logfile
    if logfile is None:
        logfile = open("../error.log", "w")
    logfile.write(message + "\n")

def clean_output():
    "Clean out old output directory"
    if os.path.isdir(output_directory):
        shutil.rmtree(output_directory)
    os.mkdir(output_directory)

def generate_test_cases(bench):
    "Generate form files for all test cases."

    begin("Generating test cases")

    # Copy form files
    if bench:
        form_directory = bench_directory
    else:
        form_directory = demo_directory
    form_files = [f for f in os.listdir(form_directory) if f.endswith(".ufl")]
    form_files.sort()
    for f in form_files:
        shutil.copy("%s/%s" % (form_directory, f), ".")
    info_green("Found %d form files" % len(form_files))

    # Generate form files for forms
    info("Generating form files for extra forms: Not implemented")

    # Generate form files for elements
    if not bench:
        from elements import elements
        info("Generating form files for extra elements (%d elements)" % len(elements))
        for (i, element) in enumerate(elements):
            open("X_Element%d.ufl" % i, "w").write("element = %s" % element)

    end()

def generate_code(args):
    "Generate code for all test cases."

    # Get a list of all files
    form_files = [f for f in os.listdir(".") if f.endswith(".ufl")]
    form_files.sort()

    begin("Generating code (%d form files found)" % len(form_files))

    # Iterate over all files
    for f in form_files:

        # Generate code
        ok = run_command("ffc %s -v -f precision=8 -fconvert_exceptions_to_warnings %s" % (" ".join(args), f))

        # Check status
        if ok:
            info_green("%s OK" % f)
        else:
            info_red("%s failed" % f)

    end()

def validate_code():
    "Validate generated code against references."

    # Get a list of all files
    header_files = [f for f in os.listdir(".") if f.endswith(".h")]
    header_files.sort()

    begin("Validating generated code (%d header files found)" % len(header_files))

    # Iterate over all files
    for f in header_files:

        # Get generated code
        generated_code = open(f).read()

        # Get reference code
        reference_file = "../references/%s" % f
        if os.path.isfile(reference_file):
            reference_code = open(reference_file).read()
        else:
            info_blue("Missing reference for %s" % f)
            continue

        # Compare with reference
        if generated_code == reference_code:
            info_green("%s OK" % f)
        else:
            info_red("%s differs" % f)
            diff = "\n".join([line for line in difflib.unified_diff(reference_code.split("\n"), generated_code.split("\n"))])
            s = "Code differs for %s, diff follows" % f
            log_error("\n" + s + "\n" + len(s)*"-")
            log_error(diff)

    end()

def build_programs(bench):
    "Build test programs for all test cases."

    # Get a list of all files
    header_files = [f for f in os.listdir(".") if f.endswith(".h")]
    header_files.sort()

    begin("Building test programs (%d header files found)" % len(header_files))

    # Set compiler options
    if bench > 0:
        info("Benchmarking activated")
        # Takes too long to build with -O2
        #compiler_options = "`pkg-config --cflags ufc-1` -Wall -Werror -O2"
        compiler_options = "`pkg-config --cflags ufc-1` -Wall -Werror"
    else:
        compiler_options = "`pkg-config --cflags ufc-1` -Wall -Werror -g"
    info("Compiler options: %s" % compiler_options)

    # Iterate over all files
    for f in header_files:

        # Generate test code
        filename = generate_test_code(f, bench)

        # Compile test code
        prefix = f.split(".h")[0]
        command = "g++ %s -o %s.bin %s.cpp" % (compiler_options, prefix, prefix)
        ok = run_command(command)

        # Check status
        if ok:
            info_green("%s OK" % prefix)
        else:
            info_red("%s failed" % prefix)

    end()

def run_programs():
    "Run generated programs."

    # Get a list of all files
    test_programs = [f for f in os.listdir(".") if f.endswith(".bin")]
    test_programs.sort()

    begin("Running generated programs (%d programs found)" % len(test_programs))

    # Iterate over all files
    for f in test_programs:

        # Compile test code
        prefix = f.split(".bin")[0]
        ok = run_command("rm -f %s.out; ./%s.bin > %s.out" % (prefix, prefix, prefix))

        # Check status
        if ok:
            info_green("%s OK" % f)
        else:
            info_red("%s failed" % f)

    end()

def validate_programs():
    "Validate generated programs against references."

    # Get a list of all files
    output_files = [f for f in os.listdir(".") if f.endswith(".out")]
    output_files.sort()

    begin("Validating generated programs (%d programs found)" % len(output_files))

    # Iterate over all files
    for f in output_files:

        # Get generated output
        generated_output = open(f).read()

        # Get reference output
        reference_file = "../references/%s" % f
        if os.path.isfile(reference_file):
            reference_output = open(reference_file).read()
        else:
            info_blue("Missing reference for %s" % f)
            continue

        # Compare with reference
        ok = True
        old = [line.split(" = ") for line in reference_output.split("\n") if " = " in line]
        new = dict([line.split(" = ") for line in generated_output.split("\n") if " = " in line])
        header = "Output differs for %s, diff follows" % f
        for (key, value) in old:

            # Check if value is present
            if not key in new:
                if ok: log_error("\n" + header + "\n" + len(header)*"-")
                log_error("%s: missing value in generated code" % key)
                ok = False
                continue

            # Extract float values
            old_values = array([float(v) for v in value.split(" ")])
            new_values = array([float(v) for v in new[key].split(" ")])

            # Check that shape is correct
            if not shape(old_values) == shape(new_values):
                if ok: log_error("\n" + header + "\n" + len(header)*"-")
                log_error("%s: shape mismatch" % key)
                ok = False
                continue

            # Check that values match to within tolerance
            diff = max(abs(old_values - new_values))
            if diff > tolerance:
                if ok: log_error("\n" + header + "\n" + len(header)*"-")
                log_error("%s: values differ, error = %g (tolerance = %g)" % (key, diff, tolerance))
                log_error("  old = " + " ".join("%.16g" % v for v in old_values))
                log_error("  new = " + " ".join("%.16g" % v for v in new_values))
                ok = False

        # Add debugging output to log file
        debug = "\n".join([line for line in generated_output.split("\n") if "debug" in line])
        if debug: log_error(debug)

        # Check status
        if ok:
            info_green("%s OK" % f)
        else:
            info_red("%s differs" % f)

    end()

def main(args):
    "Run all regression tests."

    # Clean out old output directory
    clean_output()

    # Enter output directory
    os.chdir(output_directory)

    # Check command-line arguments
    bench = "--bench" in args
    fast = "--fast" in args
    args = [arg for arg in args if not arg in ("--bench", "--fast")]

    # Generate test cases
    generate_test_cases(bench)

    # Generate code
    generate_code(args)

    # Validate code
    if not bench:
        validate_code()

    # Build, run and validate programs
    if fast:
        info("Skipping program validation")
    elif bench:
        build_programs(bench)
        run_programs()
    else:
        build_programs(bench)
        run_programs()
        validate_programs()

    # Print results
    if logfile is None:
        info_green("Regression tests OK")
        return 0
    else:
        info_red("Regression tests failed")
        info("Error messages stored in error.log")
        return 1

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
