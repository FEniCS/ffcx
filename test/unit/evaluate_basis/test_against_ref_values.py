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
# First added:  2010-01-18
# Last changed: 2010-01-29

from cppcode import evaluate_basis_code
from ufl import FiniteElement, MixedElement
from instant.output import get_status_output
import sys, os, pickle, numpy, shutil
from elements import single_elements, mixed_elements

ffc_failed = []
gcc_failed = []
run_failed = []

def check_results(values, reference):
    "Check results and print summary."

    # Check if we have missing values.
    missing_vals = []
    num_refs = len(list(reference.keys()))
    for element in list(reference.keys()):
        if not element in values:
            missing_vals.append(element)

    missing_refs = []
    diffs = []
    correct = []
    num_ok = 0
    print("")
    sorted_elements = sorted(values.keys())
    for element in sorted_elements:
        vals = values[element]
        print("\nResults for %s:" % element)

        if vals is None:
            print("Error")
            continue

        # Get reference values
        if not element in reference:
            missing_refs.append(element)
            print("Missing reference")
            continue
        refs = reference[element]
        tol = 1e-12

        e = max(abs(vals - refs))
        if e < tol:
            num_ok += 1
            print("OK: (diff = %g)" % e)
            correct.append(element)
        else:
            print("*** (diff = %g)" % e)
            diffs.append(element)

    if ffc_failed == gcc_failed == run_failed == missing_refs == diffs == missing_vals:
        print("\nAll %d elements verified OK" % len(reference))
        return 0
    else:
        print("\n*** The values were correct for the following elements:\n" + "\n\n".join(correct))
    if len(ffc_failed) > 0:
        print("\n*** FFC compilation failed for the following elements:\n" + "\n\n".join(ffc_failed))
    if len(gcc_failed) > 0:
        print("\n*** g++ compilation failed for the following elements:\n" + "\n\n".join(gcc_failed))
    if len(run_failed) > 0:
        print("\n*** Evaluation failed (seg. fault?) for the following elements:\n" + "\n\n".join(run_failed))
    if len(missing_refs) > 0:
        print("\n*** No reference values were found for the following elements:\n" + "\n\n".join(missing_refs))
    if len(missing_vals) > 0:
        print("\n*** No values were computed the following %d elements:\n" % len(missing_vals) +\
              "\n\n".join(missing_vals))
    if len(diffs) > 0:
        print("\n*** Difference in values were found for the following elements:\n" + "\n\n".join(diffs))

    num_ffc = len(ffc_failed)
    num_gcc = len(gcc_failed)
    num_run = len(run_failed)
    num_ref = len(missing_refs)
    num_val = len(missing_vals)
    num_cor = len(correct)
    num_dif = len(diffs)
    print("\nNum ref elements: ", num_refs)
    print("Num ffc fail:     ", num_ffc)
    print("Num gcc fail:     ", num_gcc)
    print("Num run fail:     ", num_run)
    print("Num miss ref:     ", num_ref)
    print("Num miss val:     ", num_val)
    print("Num correct:      ", num_cor)
    print("Num diff:         ", num_dif)
    print("Total:            ", num_ffc + num_gcc + num_run + num_ref + num_val + num_cor + num_dif)

    return 1

def compile_element(ufl_element):
    "Create UFL form file with a single element in it and compile it with FFC"
    f = open("test.ufl", "w")
    if isinstance(ufl_element, (FiniteElement, MixedElement)):
        f.write("element = " + repr(ufl_element))
    f.close()
    error, out = get_status_output("ffc test.ufl")
    if error:
        ffc_failed.append(repr(ufl_element))
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

def compute_values(ufl_element):
    "Compute values of basis functions for given element."

    # Get relevant element name
    element_name = get_element_name(ufl_element)

    # Create g++ code
    options = {"element": element_name}
    code = evaluate_basis_code % options
    f = open("evaluate_basis.cpp", "w")
    f.write(code)
    f.close()

    # Get UFC flags
    ufc_cflags = get_status_output("pkg-config --cflags ufc-1")[1].strip()

    # Compile g++ code
    c = "g++ %s -Wall -Werror -o evaluate_basis evaluate_basis.cpp" % ufc_cflags
    error, output = get_status_output(c)
    if error:
        gcc_failed.append(repr(ufl_element))
        return None

    # Run compiled code and get values
    error, output = get_status_output(".%sevaluate_basis" % os.path.sep)
    if error:
        run_failed.append(repr(ufl_element))
        return None
    values = [float(value) for value in output.split(" ") if len(value) > 0]
    return numpy.array(values)

def print_refs():
    if os.path.isfile("reference.pickle"):
        reference = pickle.load(open("reference.pickle", "r"))
        for elem, vals in list(reference.items()):
            print()
            print(elem)
            print(vals)
    else:
        raise RuntimeError("No references to print")

def main(args):
    "Call evaluate basis for a range of different elements."

    if "refs" in args:
        print_refs()
        return 0

    # Change to temporary folder and copy form files
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")
    os.chdir("tmp")

    values = {}
    # Evaluate basis for single elements
    print("\nComputing evaluate_basis for single elements")
    for element in single_elements:
        for shape in element["shapes"]:
            for order in element["orders"]:
                ufl_element = FiniteElement(element["family"], shape, order)
                print("Compiling element: ", str(ufl_element))
                error = compile_element(ufl_element)
                if error:
                    values[repr(ufl_element)] = None
                    continue
                print("Computing values")
                values[repr(ufl_element)] = compute_values(ufl_element)

    # Evaluate basis for mixed elements
    print("\nComputing evaluate_basis for mixed elements")
    for ufl_element in mixed_elements:
        print("Compiling element: ", str(ufl_element))
        error = compile_element(ufl_element)
        if error:
            values[repr(ufl_element)] = None
            continue
        print("Computing values")
        values[repr(ufl_element)] = compute_values(ufl_element)

    # Load or update reference values
    os.chdir(os.pardir)
    if os.path.isfile("reference.pickle"):
        reference = pickle.load(open("reference.pickle", "r"))
    else:
        print("Unable to find reference values, storing current values.")
        pickle.dump(values, open("reference.pickle", "w"))
        return 0

    # Check results
    error = check_results(values, reference)

    if not error:
        # Remove temporary directory
        shutil.rmtree("tmp")

    return error

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
