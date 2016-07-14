# Copyright (C) 2010 Kristian B. Oelgaard
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
# First added:  2010-01-29
# Last changed: 2013-02-06

from ffc.log import info, info_red, info_blue, info_green, debug
from instant.output import get_status_output
import numpy
import os
import sys
from six.moves import xrange as range

tol = 1e-14
crit_tol = 1e-8

def xcomb(items, n):
    "Create n-tuples with combinations of items."
    if n==0: yield []
    else:
        for i in range(len(items)):
            for cc in xcomb(items[:i]+items[i+1:],n-1):
                yield [items[i]]+cc


# Global log file
def log_error(message, log_file):
    "Log error message."
    log = open(os.path.join(os.pardir, log_file), "a")
    log.write("\n" + "-"*79 + "\n" + message + "\n" + "-"*79 + "\n")
    log.close()


def print_results(num_tests, ffc_fail, gcc_fail, run_fail, dif_cri, dif_acc,
                  correct):
    "Check print summary."

    num_ffc = len(ffc_fail)
    num_gcc = len(gcc_fail)
    num_run = len(run_fail)
    num_cri = len(dif_cri)
    num_acc = len(dif_acc)
    num_cor = len(correct)

    if ffc_fail == gcc_fail == run_fail == dif_cri == dif_acc == []:
        info_green("\nAll %d tests OK" % num_tests)
        return 0

    num_tests = str(num_tests)
    num_tot = str(num_ffc + num_gcc + num_run + num_cor + num_cri + num_acc)
    num_ffc = str(num_ffc)
    num_gcc = str(num_gcc)
    num_run = str(num_run)
    num_cor = str(num_cor)
    num_cri = str(num_cri)
    num_acc = str(num_acc)
    num_ffc = " "*(len(num_tests) - len(num_ffc)) + num_ffc
    num_gcc = " "*(len(num_tests) - len(num_gcc)) + num_gcc
    num_run = " "*(len(num_tests) - len(num_run)) + num_run
    num_cor = " "*(len(num_tests) - len(num_cor)) + num_cor
    num_cri = " "*(len(num_tests) - len(num_cri)) + num_cri
    num_acc = " "*(len(num_tests) - len(num_acc)) + num_acc
    num_tot = " "*(len(num_tests) - len(num_tot)) + num_tot

    info("\n\n*************** SUMMARY ***************")
    info("\n  Number of tests:                  " + num_tests)
    info("\n  Num ffc fail:                     " + num_ffc)
    info("  Num gcc fail:                     " + num_gcc)
    info("  Num run fail:                     " + num_run)
    info(("  Num correct:        (tol. %g): " % tol) + num_cor)
    info(("  Num diff. critical: (tol. %g): " % crit_tol) + num_cri)
    info("  Num diff. acceptable:             " + num_acc)
    info("  Total:                            " + num_tot)
    info("")
    # Return 0 if there was only acceptable errors.
    if ffc_fail == gcc_fail == run_fail == dif_cri == []:
        return 0
    return 1


def compile_element(ufl_element, ffc_fail, log_file):
    "Create UFL form file with a single element in it and compile it with FFC"
    f = open("test.ufl", "w")
    f.write("element = " + repr(ufl_element))
    f.close()
    error, output = get_status_output("ffc test.ufl")
    if error:
        info_red("FFC compilation failed.")
        log_error("element: %s,\n%s\n" % (str(ufl_element), output), log_file)
        ffc_fail.append(str(ufl_element))
    return error


def get_element_name(ufl_element):
    "Extract relevant element name from header file."
    f = open("test.h")
    lines = f.readlines()
    f.close()

    signature = repr(ufl_element)
    name = None
    for e, l in enumerate(lines):
        if "class" in l and "finite_element" in l:
            name = l
        if signature in l:
            break
    if name is None:
        raise RuntimeError("No finite element class found")
    return name.split()[1][:-1]


def compile_gcc_code(ufl_element, code, gcc_fail, log_file):
    # Write code.
    f = open("evaluate_basis.cpp", "w")
    f.write(code)
    f.close()

    # Get location of UFC file
    import ffc.ufc_include
    ufc_include_path = ffc.ufc_include.get_ufc_include()

    # Compile g++ code
    c = "g++ -I{} -std=c++11 -Wall -Werror -o evaluate_basis_test_code evaluate_basis.cpp".format(ufc_include_path)
    f = open("compile.sh", "w")
    f.write(c + "\n")
    f.close()
    error, output = get_status_output(c)
    if error:
        info_red("GCC compilation failed.")
        log_error("element: %s,\n%s\n" % (str(ufl_element), output), log_file)
        gcc_fail.append(str(ufl_element))
        if error and ("-f" in sys.argv or "--failfast" in sys.argv):
            print("FAIL")
            exit(1)
        return error
    else:
        return None


def run_code(ufl_element, deriv_order, run_fail, log_file):
    "Compute values of basis functions for given element."

    # Run compiled code and get values
    error, output = get_status_output(".%sevaluate_basis_test_code %d" % (os.path.sep,
                                                                deriv_order))
    if error:
        info_red("Runtime error (segmentation fault?).")
        log_error("element: %s,\n%s\n" % (str(ufl_element), output), log_file)
        run_fail.append(str(ufl_element))
        return None
    values = [[float(value) for value in line.strip().split(" ") if value] for line in output.strip().split("\n")]
    return numpy.array(values)

def verify_values(ufl_element, ref_values, ffc_values, dif_cri, dif_acc, correct, log_file):
    "Check the values from evaluate_basis*() against some reference values."

    num_tests = len(ffc_values)
    if num_tests != len(ref_values):
        raise RuntimeError("The number of computed values is not equal to the number of reference values.")

    errors = [str(ufl_element)]
    for deriv_order in range(num_tests):
        s = ""
        if deriv_order == 0:
            s = "  evaluate_basis"
        else:
            s = "  evaluate_basis_derivatives, order = %d" % deriv_order
        e = abs(ffc_values[deriv_order] - ref_values[deriv_order])
        error = e.max()
        if error > tol:
            if error >  crit_tol:
                m = "%s failed: error = %s (crit_tol: %s)" % (s, str(error), str(crit_tol))
                info_red(m)
                dif_cri.append(str(ufl_element))
                s = s + "\n" + m
            else:
                m = "%s ok: error = %s (tol: %s)" % (s, str(error), str(tol))
                info_blue(m)
                dif_acc.append(str(ufl_element))
                s = s + "\n" + m
            errors.append(s)
        else:
            info_green("%s OK" % s)
            correct.append(str(ufl_element))
    # Log errors if any
    if len(errors) > 1:
        log_error("\n".join(errors), log_file)

    return num_tests
