# -*- coding: utf-8 -*-
"""This script compiles and verifies the output for all form files
found in the 'demo' directory. The verification is performed in two
steps. First, the generated code is compared with stored references.
Then, the output from all functions in the generated code is compared
with stored reference values.

This script can also be used for benchmarking tabulate_tensor for all
form files found in the 'bench' directory. To run benchmarks, use the
option --bench.

"""

# Copyright (C) 2010-2013 Anders Logg, Kristian B. Oelgaard and Marie E. Rognes
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
# Modified by Martin Sandve Alnæs, 2013-2017
# Modified by Johannes Ring, 2013
# Modified by Kristian B. Oelgaard, 2013
# Modified by Garth N. Wells, 2014

# FIXME: Need to add many more test cases. Quite a few DOLFIN forms
# failed after the FFC tests passed.

from collections import OrderedDict, defaultdict
import os
import sys
import shutil
import difflib
import sysconfig
import subprocess
import time
import logging
import traceback
from numpy import array, shape, abs, max, isnan
import ffc
from ffc.log import begin, end, info, info_red, info_green, info_blue
from ffc.log import ffc_logger, ERROR
from ufl.log import ufl_logger
from ufl.utils.str import as_native_str
from ffc import get_ufc_cxx_flags
from ffc.backends.ufc import get_include_path as get_ufc_include
from ufctest import generate_test_code

# Parameters
output_tolerance = 1e-5
demo_directory = "../../../../demo"
bench_directory = "../../../../bench"

# Global log file
logfile = "error.log"

# Remove old log file
if os.path.isfile(logfile):
    os.remove(logfile)

class GEFilter(object):
    """Filter messages that are greater or equal to given log level"""
    def __init__(self, level):
        self.__level = level

    def filter(self, record):
        return record.levelno >= self.__level

class LTFilter(object):
    """Filter messages that are less than given log level"""
    def __init__(self, level):
        self.__level = level

    def filter(self, record):
        return record.levelno < self.__level


# Filter out error messages from std output
splitlevel = ERROR
ffc_logger.get_handler().addFilter(LTFilter(splitlevel))
ufl_logger.get_handler().addFilter(LTFilter(splitlevel))

# Filter out error messages to log file
file_handler = logging.FileHandler(logfile)
file_handler.addFilter(GEFilter(splitlevel))
ffc_logger.get_logger().addHandler(file_handler)
ufl_logger.get_logger().addHandler(file_handler)

# Extended quadrature tests (optimisations)
ext_quad = [
    "-r quadrature -O -feliminate_zeros",
    "-r quadrature -O -fsimplify_expressions",
    "-r quadrature -O -fprecompute_ip_const",
    "-r quadrature -O -fprecompute_basis_const",
    "-r quadrature -O -fprecompute_ip_const -feliminate_zeros",
    "-r quadrature -O -fprecompute_basis_const -feliminate_zeros",
]

# Extended uflacs tests
# (to be extended with optimisation parameters later)
ext_uflacs = [
    "-r uflacs -O -fvectorize -fpadlen=4 -falignas=32",
    "-r uflacs -O -fno-enable_sum_factorization",
    "-r uflacs -O -fno-enable_preintegration",
    "-r uflacs -O -fenable_premultiplication",
]

known_quad_failures = set([
    "PoissonQuad.ufl",
])

known_uflacs_failures = set([
    "CustomIntegral.ufl",
    "CustomMixedIntegral.ufl",
    "CustomVectorIntegral.ufl",
    "MetaData.ufl",
])

known_tsfc_failures = set([
    # Expected not to work
    "CustomIntegral.ufl",
    "CustomMixedIntegral.ufl",
    "CustomVectorIntegral.ufl",
    "MetaData.ufl",
])


_command_timings = []


def run_command(command, verbose):
    "Run command and collect errors in log file."
    global _command_timings

    t1 = time.time()
    try:
        output = as_native_str(subprocess.check_output(command, shell=True))
        t2 = time.time()
        _command_timings.append((command, t2 - t1))
        if verbose:
            print(output)
        return True
    except subprocess.CalledProcessError as e:
        t2 = time.time()
        _command_timings.append((command, t2 - t1))
        if e.output:
            log_error(e.output)
            print(e.output)
        return False


def log_error(message):
    "Log error message."
    ffc_logger.get_logger().error(message)


def clean_output(output_directory):
    "Clean out old output directory"
    if os.path.isdir(output_directory):
        shutil.rmtree(output_directory)
    os.mkdir(output_directory)


def generate_test_cases(bench, only_forms, skip_forms):
    "Generate form files for all test cases."

    begin("Generating test cases")

    # Copy form files
    if bench:
        form_directory = bench_directory
    else:
        form_directory = demo_directory

    # Make list of form files
    form_files = [f for f in os.listdir(form_directory)
                  if f.endswith(".ufl") and f not in skip_forms]
    if only_forms:
        form_files = [f for f in form_files if f in only_forms]
    form_files.sort()

    for f in form_files:
        shutil.copy(os.path.join(form_directory, f), ".")
    info_green("Found %d form files" % len(form_files))

    # Generate form files for forms
    #info("Generating form files for extra forms: Not implemented")

    # Generate form files for elements
    if not (bench or only_forms):
        from elements import elements
        info("Generating form files for extra elements (%d elements)"
             % len(elements))
        for (i, element) in enumerate(elements):
            with open("X_Element%d.ufl" % i, "w") as f:
                f.write("element = %s" % element)

    end()


def generate_code(args, only_forms, skip_forms, debug):
    "Generate code for all test cases."
    global _command_timings

    # Get a list of all files
    form_files = [f for f in os.listdir(".")
                  if f.endswith(".ufl") and f not in skip_forms]
    if only_forms:
        form_files = [f for f in form_files if f in only_forms]
    form_files.sort()

    begin("Generating code (%d form files found)" % len(form_files))

    # TODO: Parse additional options from .ufl file? I.e. grep for
    # some sort of tag like '#ffc: <flags>'.
    special = {"AdaptivePoisson.ufl": "-e", }

    failures = []

    # Iterate over all files
    for f in form_files:
        options = [special.get(f, "")]
        options.extend(args)
        options.extend(["-f", "precision=8", "-f", "epsilon=1e-7",  "-fconvert_exceptions_to_warnings"])
        options.append(f)
        options = list(filter(None, options))

        cmd = sys.executable + " -m ffc " + " ".join(options)

        # Generate code
        t1 = time.time()
        try:
            ok = ffc.main(options)
        except Exception as e:
            if debug:
                raise e
            msg = traceback.format_exc()
            log_error(cmd)
            log_error(msg)
            ok = 1
        finally:
            t2 = time.time()
            _command_timings.append((cmd, t2 - t1))

        # Check status
        if ok == 0:
            info_green("%s OK" % f)
        else:
            info_red("%s failed" % f)
            failures.append(f)

    end()
    return failures


def validate_code(reference_dir):
    "Validate generated code against references."

    # Get a list of all files
    header_files = sorted([f for f in os.listdir(".") if f.endswith(".h")])

    begin("Validating generated code (%d header files found)"
          % len(header_files))

    failures = []

    # Iterate over all files
    for f in header_files:

        # Get generated code
        generated_code = open(f).read()

        # Get reference code
        reference_file = os.path.join(reference_dir, f)
        if os.path.isfile(reference_file):
            reference_code = open(reference_file).read()
        else:
            info_blue("Missing reference for %s" % reference_file)
            continue

        # Compare with reference
        if generated_code == reference_code:
            info_green("%s OK" % f)
        else:
            info_red("%s differs" % f)
            difflines = difflib.unified_diff(
                reference_code.split("\n"),
                generated_code.split("\n"))
            diff = "\n".join(difflines)
            s = ("Code differs for %s, diff follows (reference first, generated second)"
                 % os.path.join(*reference_file.split(os.path.sep)[-3:]))
            log_error("\n" + s + "\n" + len(s) * "-")
            log_error(diff)
            failures.append(f)

    end()
    return failures


def find_boost_cflags():
    # Get Boost dir (code copied from ufc/src/utils/python/ufc_utils/build.py)
    # Set a default directory for the boost installation
    if sys.platform == "darwin":
        # Use Brew as default
        default = os.path.join(os.path.sep, "usr", "local")
    else:
        default = os.path.join(os.path.sep, "usr")

    # If BOOST_DIR is not set use default directory
    boost_inc_dir = ""
    boost_lib_dir = ""
    boost_math_tr1_lib = "boost_math_tr1"
    boost_dir = os.getenv("BOOST_DIR", default)
    boost_is_found = False
    for inc_dir in ["", "include"]:
        if os.path.isfile(os.path.join(boost_dir, inc_dir, "boost",
                                       "version.hpp")):
            boost_inc_dir = os.path.join(boost_dir, inc_dir)
            break
    libdir_multiarch = "lib/" + sysconfig.get_config_vars().get("MULTIARCH", "")
    for lib_dir in ["", "lib", libdir_multiarch, "lib64"]:
        for ext in [".so", "-mt.so", ".dylib", "-mt.dylib"]:
            _lib = os.path.join(boost_dir, lib_dir, "lib" + boost_math_tr1_lib
                                + ext)
            if os.path.isfile(_lib):
                if "-mt" in _lib:
                    boost_math_tr1_lib += "-mt"
                boost_lib_dir = os.path.join(boost_dir, lib_dir)
                break
    if boost_inc_dir != "" and boost_lib_dir != "":
        boost_is_found = True

    if boost_is_found:
        boost_cflags = " -I%s -L%s" % (boost_inc_dir, boost_lib_dir)
        boost_linkflags = "-l%s" % boost_math_tr1_lib
    else:
        boost_cflags = ""
        boost_linkflags = ""
        info_red("""The Boost library was not found.
If Boost is installed in a nonstandard location,
set the environment variable BOOST_DIR.
Forms using bessel functions will fail to build.
""")
    return boost_cflags, boost_linkflags


def build_programs(bench, permissive, debug, verbose):
    "Build test programs for all test cases."

    # Get a list of all files
    header_files = sorted([f for f in os.listdir(".") if f.endswith(".h")])

    begin("Building test programs (%d header files found)" % len(header_files))

    # Get UFC flags
    ufc_cflags = "-I" + get_ufc_include() + " " + " ".join(get_ufc_cxx_flags())

    # Get boost flags
    boost_cflags, boost_linkflags = find_boost_cflags()

    # Get compiler
    compiler = os.getenv("CXX", "g++")

    # Set compiler options
    compiler_options = " -Wall"
    if not permissive:
        compiler_options += " -Werror -pedantic"

    # Always need ufc
    compiler_options += " " + ufc_cflags

    if bench:
        info("Benchmarking activated")
        compiler_options += " -O3 -march=native"
        # Workaround for gcc bug: gcc is too eager to report array-bounds warning with -O3
        compiler_options += " -Wno-array-bounds"

    if debug:
        info("Debugging activated")
        compiler_options += " -g -O0"

    info("Compiler options: %s" % compiler_options)

    failures = []

    # Iterate over all files
    for f in header_files:
        prefix = f.split(".h")[0]

        # Options for all files
        cpp_flags = compiler_options
        ld_flags = ""

        # Only add boost flags if necessary
        needs_boost = prefix == "MathFunctions"
        if needs_boost:
            info("Additional compiler options for %s: %s" % (prefix, boost_cflags))
            info("Additional linker options for %s: %s" % (prefix, boost_linkflags))
            cpp_flags += " " + boost_cflags
            ld_flags += " " + boost_linkflags

        # Generate test code
        filename = generate_test_code(f)

        # Compile test code
        command = "%s %s -o %s.bin %s.cpp %s" % \
                  (compiler, cpp_flags, prefix, prefix, ld_flags)
        ok = run_command(command, verbose)

        # Store compile command for easy reproduction
        with open("%s.build" % (prefix,), "w") as f:
            f.write(command + "\n")

        # Check status
        if ok:
            info_green("%s OK" % prefix)
        else:
            info_red("%s failed" % prefix)
            failures.append(prefix)

    end()
    return failures


def run_programs(bench, debug, verbose):
    "Run generated programs."

    # This matches argument parsing in the generated main files
    bench = 'b' if bench else ''

    # Get a list of all files
    test_programs = sorted([f for f in os.listdir(".") if f.endswith(".bin")])

    begin("Running generated programs (%d programs found)" % len(test_programs))

    failures = []

    # Iterate over all files
    for f in test_programs:
        # Compile test code
        prefix = f.split(".bin")[0]
        ok = run_command(".%s%s.bin %s" % (os.path.sep, prefix, bench), verbose)

        # Check status
        if ok:
            info_green("%s OK" % f)
        else:
            info_red("%s failed" % f)
            failures.append(f)
    end()
    return failures


def validate_programs(reference_dir):
    "Validate generated programs against references."

    # Get a list of all files
    output_files = sorted(f for f in os.listdir(".") if f.endswith(".json"))

    begin("Validating generated programs (%d .json program output files found)"
          % len(output_files))

    failures = []

    # Iterate over all files
    for fj in output_files:

        # Get generated json output
        if os.path.exists(fj):
            generated_json_output = open(fj).read()
            if "nan" in generated_json_output:
                info_red("Found nan in generated json output, replacing with 999 to be able to parse as python dict.")
                generated_json_output = generated_json_output.replace("nan",
                                                                      "999")
        else:
            generated_json_output = "{}"

        # Get reference json output
        reference_json_file = os.path.join(reference_dir, fj)
        if os.path.isfile(reference_json_file):
            reference_json_output = open(reference_json_file).read()
        else:
            info_blue("Missing reference for %s" % reference_json_file)
            reference_json_output = "{}"

        # Compare json with reference using recursive diff algorithm
        # TODO: Write to different error file?
        from recdiff import recdiff, print_recdiff, DiffEqual
        # Assuming reference is well formed
        reference_json_output = eval(reference_json_output)
        try:
            generated_json_output = eval(generated_json_output)
        except Exception as e:
            info_red("Failed to evaluate json output for %s" % fj)
            log_error(str(e))
            generated_json_output = None
        json_diff = (None if generated_json_output is None else
                     recdiff(generated_json_output, reference_json_output, tolerance=output_tolerance))
        json_ok = json_diff == DiffEqual

        # Check status
        if json_ok:
            info_green("%s OK" % fj)
        else:
            info_red("%s differs" % fj)
            log_error("Json output differs for %s, diff follows (generated first, reference second)"
                      % os.path.join(*reference_json_file.split(os.path.sep)[-3:]))
            print_recdiff(json_diff, printer=log_error)
            failures.append(fj)

    end()
    return failures


def main(args):
    "Run all regression tests."

    # Check command-line arguments TODO: Use argparse
    only_auto  = "--only-auto" in args
    use_auto   = "--skip-auto" not in args
    use_uflacs = "--skip-uflacs" not in args
    use_quad   = "--skip-quad" not in args
    use_tsfc   = "--use-tsfc" in args
    use_ext_quad   = "--ext-quad" in args
    use_ext_uflacs = "--ext-uflacs" in args

    skip_download = "--skip-download" in args
    skip_run = "--skip-run" in args
    skip_code_diff = "--skip-code-diff" in args
    skip_validate = "--skip-validate" in args
    bench = "--bench" in args
    debug = "--debug" in args
    verbose = ("--verbose" in args) or debug  # debug implies verbose

    permissive = "--permissive" in args or bench
    tolerant = "--tolerant" in args
    print_timing = "--print-timing" in args
    show_help = "--help" in args

    flags = (
        "--only-auto",
        "--skip-auto",
        "--skip-uflacs",
        "--skip-quad",
        "--use-tsfc",
        "--ext-quad",
        "--skip-download",
        "--skip-run",
        "--skip-code-diff",
        "--skip-validate",
        "--bench",
        "--debug",
        "--verbose",
        "--permissive",
        "--tolerant",
        "--print-timing",
        "--help",
    )
    args = [arg for arg in args if arg not in flags]

    # Hack: add back --verbose for ffc.main to see
    if verbose:
        args = args + ["--verbose"]

    if show_help:
        info("Valid arguments:\n" + "\n".join(flags))
        return 0

    if bench or not skip_validate:
        skip_run = False
    if bench:
        skip_code_diff = True
        skip_validate = True
    if use_ext_quad or use_ext_uflacs:
        skip_code_diff = True

    # Extract .ufl names from args
    only_forms = set([arg for arg in args if arg.endswith(".ufl")])
    args = [arg for arg in args if arg not in only_forms]

    # Download reference data
    if skip_download:
        info_blue("Skipping reference data download")
    else:
        try:
            cmd = "./scripts/download"
            output = as_native_str(subprocess.check_output(cmd, shell=True))
            print(output)
            info_green("Download reference data ok")
        except subprocess.CalledProcessError as e:
            print(e.output)
            info_red("Download reference data failed")

    if tolerant:
        global output_tolerance
        output_tolerance = 1e-3

    # Clean out old output directory
    output_directory = "output"
    clean_output(output_directory)
    os.chdir(output_directory)

    # Adjust which test cases (combinations of compile arguments) to run here
    test_cases = []
    if only_auto:
        test_cases += ["-r auto"]
    else:
        if use_auto:
            test_cases += ["-r auto"]
        if use_uflacs:
            test_cases += ["-r uflacs -O0", "-r uflacs -O"]
        if use_quad:
            test_cases += ["-r quadrature -O0", "-r quadrature -O"]
            import warnings
            from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning
            warnings.simplefilter("once", QuadratureRepresentationDeprecationWarning)
        if use_tsfc:
            test_cases += ["-r tsfc -O0", "-r tsfc -O"]
            # Silence good-performance messages by COFFEE
            import coffee
            coffee.set_log_level(coffee.logger.PERF_WARN)
        if use_ext_quad:
            test_cases += ext_quad
        if use_ext_uflacs:
            test_cases += ext_uflacs

    test_case_timings = {}

    fails = OrderedDict()

    for argument in test_cases:
        test_case_timings[argument] = time.time()
        fails[argument] = OrderedDict()

        begin("Running regression tests with %s" % argument)

        # Clear and enter output sub-directory
        sub_directory = "_".join(argument.split(" ")).replace("-", "")
        clean_output(sub_directory)
        os.chdir(sub_directory)

        # Workarounds for feature lack in representation
        if "quadrature" in argument and not only_forms:
            skip_forms = known_quad_failures
            info_blue("Skipping forms known to fail with quadrature:\n" + "\n".join(sorted(skip_forms)))
        elif "uflacs" in argument and not only_forms:
            skip_forms = known_uflacs_failures
            info_blue("Skipping forms known to fail with uflacs:\n" + "\n".join(sorted(skip_forms)))
        elif "tsfc" in argument and not only_forms:
            skip_forms = known_tsfc_failures
            info_blue("Skipping forms known to fail with tsfc:\n" + "\n".join(sorted(skip_forms)))
        else:
            skip_forms = set()

        # Generate test cases
        generate_test_cases(bench, only_forms, skip_forms)

        # Generate code
        failures = generate_code(args + argument.split(), only_forms, skip_forms, debug)
        if failures:
            fails[argument]["generate_code"] = failures

        # Location of reference directories
        reference_directory = os.path.abspath("../../ffc-reference-data/")
        code_reference_dir = os.path.join(reference_directory, sub_directory)

        # Note: We use the r_auto references for all test cases. This
        # ensures that we continously test that the codes generated by
        # all different representations are equivalent.
        output_reference_dir = os.path.join(reference_directory, "r_auto")

        # Validate code by comparing to code generated with this set
        # of compiler parameters
        if skip_code_diff:
            info_blue("Skipping code diff validation")
        else:
            failures = validate_code(code_reference_dir)
            if failures:
                fails[argument]["validate_code"] = failures

        # Build and run programs and validate output to common
        # reference
        if skip_run:
            info_blue("Skipping program execution")
        else:
            failures = build_programs(bench, permissive, debug, verbose)
            if failures:
                fails[argument]["build_programs"] = failures

            failures = run_programs(bench, debug, verbose)
            if failures:
                fails[argument]["run_programs"] = failures

            # Validate output to common reference results
            if skip_validate:
                info_blue("Skipping program output validation")
            else:
                failures = validate_programs(output_reference_dir)
                if failures:
                    fails[argument]["validate_programs"] = failures

        # Go back up
        os.chdir(os.path.pardir)

        end()
        test_case_timings[argument] = time.time() - test_case_timings[argument]

    # Go back up
    os.chdir(os.path.pardir)

    # Print results
    if print_timing:
        info_green("Timing of all commands executed:")
        timings = '\n'.join("%10.2e s  %s" % (t, name) for (name, t)
                            in _command_timings)
        info_blue(timings)

    for argument in test_cases:
        info_blue("Total time for %s: %.1f s" % (argument, test_case_timings[argument]))

    num_failures = sum(len(failures_phase)
                       for failures_args in fails.values()
                       for failures_phase in failures_args.values())

    if num_failures == 0:
        info_green("Regression tests OK")
        return 0
    else:
        info_red("Regression tests failed")
        info_red("")
        info_red("Long summary:")
        for argument in test_cases:
            if not fails[argument]:
                info_green("  No failures with args '%s'" % argument)
            else:
                info_red("  Failures with args '%s':" % argument)
                for phase, failures in fails[argument].items():
                    info_red("    %d failures in %s:" % (len(failures), phase))
                    for f in failures:
                        info_red("      %s" % (f,))
        info_red("")
        info_red("Short summary:")
        phase_fails = defaultdict(int)
        for argument in test_cases:
            if not fails[argument]:
                info_green("  No failures with args '%s'" % argument)
            else:
                info_red("  Number of failures with args '%s':" % argument)
                for phase, failures in fails[argument].items():
                    info_red("    %d failures in %s." % (len(failures), phase))
                    phase_fails[phase] += len(failures)
        info_red("")
        info_red("Total failures for all args:")
        for phase, count in phase_fails.items():
            info_red("    %s: %d failed" % (phase, count))
        info_red("")
        info_red("Error messages stored in %s" % logfile)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
