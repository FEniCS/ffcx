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
# Last changed: 2013-01-31

from .cppcode import evaluate_basis_code_fiat
from ufl import FiniteElement
from ffc.fiatinterface import create_element, reference_cell
from ffc.mixedelement import MixedElement as FFCMixedElement
from ffc.log import info, error, debug
import numpy
import sys, os, numpy, shutil
from .test_common import compile_element, print_results, compile_gcc_code,\
                        run_code, get_element_name, verify_values
import time
from ffc.log import push_level, pop_level, CRITICAL, INFO
from .elements import single_elements, mixed_elements

# Some random points
random_points = {1: [(0.114,), (0.349,), (0.986,)],
                  2: [(0.114, 0.854), (0.349, 0.247), (0.986, 0.045)],
                  3: [(0.114, 0.854, 0.126), (0.349, 0.247, 0.457), (0.986, 0.045, 0.127)]}

ffc_fail = []
gcc_fail = []
run_fail = []
dif_cri = []
dif_acc = []
correct = []
log_file = "fiat_errors.log"
def matrix(points):
    return "{%s};" % ", ".join(["{%s}" % ", ".join([str(c) for c in p]) for p in points])

def get_data(ufl_element):
    "Get needed data to run tests."

    # Create fiat element.
    element = create_element(ufl_element)

    # The derivative order that we are interested in is the degree of the element.
    if isinstance(element, FFCMixedElement):
        deriv_order = max([e.degree() for e in element.elements()])
    else:
        deriv_order = element.degree()

    # Get coordinates of the reference cell.
    cell = ufl_element.cell()
    ref_coords = reference_cell(cell.cellname()).get_vertices()

    # Get the locations of the fiat element dofs.
    elem_points =  [list(L.pt_dict.keys())[0] for L in element.dual_basis()]

    # Add some random points.
    geo_dim = cell.geometric_dimension()
    points = elem_points + random_points[geo_dim]

    return (element, points, geo_dim, ref_coords, deriv_order)

def compute_derivative_combinations(deriv_order, geo_dim):
    "Compute combinations of derivatives in spatial directions (like code snippet)."

    if deriv_order == 0:
        return [(0,)*geo_dim]

    num_derivatives = geo_dim**deriv_order
    combinations = [[0]*deriv_order for n in range(num_derivatives)]

    for i in range(1, num_derivatives):
        for k in range(i):
            j = deriv_order - 1
            while j + 1 > 0:
                j -= 1
                if combinations[i][j] + 1 > geo_dim - 1:
                    combinations[i][j] = 0
                else:
                    combinations[i][j] += 1
                    break
    # Convert to fiat tuples.
    for i in range(num_derivatives):
        combinations[i] = to_fiat_tuple(combinations[i], geo_dim)

    return combinations

def to_fiat_tuple(comb, geo_dim):
    """Convert a list of combinations of derivatives to a fiat tuple of derivatives.
    FIAT expects a list with the number of derivatives in each spatial direction.
    E.g., in 2D: u_{xyy} --> [0, 1, 1] in FFC --> (1, 2) in FIAT."""
    new_comb = [0]*geo_dim
    if comb == []:
        return tuple(new_comb)
    for i in range(geo_dim):
        new_comb[i] = comb.count(i)
    return tuple(new_comb)

def get_fiat_values(ufl_element):
    """Create a FIAT element and use it to tabulate the values on the reference
    element. The return values is a dictionary with keys equal to the derivative
    order and values is a matrix where each row is the basis values at a point.
    E.g., {0:[[1,0,0],[0,1,0], [0,0,1]]}."""

    # Get data and tabulate basis values.
    element, points, geo_dim, ref_coords, deriv_order = get_data(ufl_element)
    values = element.tabulate(deriv_order, points)
    return_values = {}
    value_shape = element.value_shape()

    # Rearrange values to match what we get from evaluate_basis*()
    for n in range(deriv_order + 1):
        combinations = compute_derivative_combinations(n, geo_dim)
        vals = []
        # If we're evaluating the basis functions, use all points, but if we're
        # evaluating the derivatives, just use the 3 arbitrary points to avoid
        # the number of tests exploding with spacedim**2.
        if n == 0:
            new_points = points
        else:
            new_points = points[-3:]
        for p, point in enumerate(new_points):
            if n != 0:
                p += element.space_dimension()
            row = [[] for i in range(element.space_dimension())]
            for i in range(element.space_dimension()):
                if value_shape == ():
                    for deriv in combinations:
                        deriv_vals = values[deriv]
                        row[i].append(deriv_vals[i][p])
                elif len(value_shape) == 1:
                    for c in range(element.value_shape()[0]):
                        for deriv in combinations:
                            deriv_vals = values[deriv]
                            row[i].append(deriv_vals[i][c][p])
                elif len(value_shape) == 2:
                    for j in range(element.value_shape()[0]):
                        for k in range(element.value_shape()[1]):
                            for deriv in combinations:
                                deriv_vals = values[deriv]
                                row[i].append(deriv_vals[i][j][k][p])
                else:
                    print(values)
                    error("Did not expect tensor elements of rank > 2")
            new_row = []
            for r in row:
                new_row += r
            vals.append(new_row)
        return_values[n] = numpy.array(vals)
    return return_values

def get_ffc_values(ufl_element):
    "Get the values from evaluate_basis and evaluate_basis_derivatives."

    # Get data and tabulate basis values.
    element, points, geo_dim, ref_coords, deriv_order = get_data(ufl_element)

    # Get relevant element name.
    element_name = get_element_name(ufl_element)

    # Create g++ code and compile.
    num_coords = len(ref_coords)
    options = {"element": element_name,
               "dim": geo_dim,
               "num_points": len(points),
               "points": matrix(points),
               "cell_ref_coords": "double cell_ref_coords[%d][%d] = %s" % (num_coords, geo_dim, matrix(ref_coords)),
               "num_coords": num_coords}
    error = compile_gcc_code(ufl_element, evaluate_basis_code_fiat % options, gcc_fail, log_file)
    if error:
        return None

    # Loop derivative order and compute values.
    ffc_values = {}
    for n in range(deriv_order + 1):
        values = run_code(ufl_element, n, run_fail, log_file)
        if values is None:
            return None
        ffc_values[n] = values
    return ffc_values

def verify_element(num_elements, i, ufl_element):
    info("\nVerifying element %d of %d: %s" % (i, num_elements, str(ufl_element)))
    error = compile_element(ufl_element, ffc_fail, log_file)

    # Return if test failed
    if error:
        return 1

    # Get FIAT values that are formatted in the same way as the values from
    # evaluate_basis and evaluate_basis_derivatives.
    # t = time.time()
    fiat_values = get_fiat_values(ufl_element)
    # print "fiat_vals: ", time.time() - t

    # Get FFC values.
    t = time.time()
    ffc_values = get_ffc_values(ufl_element)
    if ffc_values is None:
        return 1
    debug("  time to compute FFC values: %f" % (time.time() - t))

    # Compare values and return number of tests.
    return verify_values(ufl_element, fiat_values, ffc_values, dif_cri, dif_acc, correct, log_file)

def main(debug_level):
    "Call evaluate basis for a range of different elements."

    push_level(debug_level)

    # Remove old log file.
    if os.path.isfile(log_file):
        os.remove(log_file)

    # Change to temporary folder and copy form files
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")
    os.chdir("tmp")

    # Create list of all elements that have to be tested.
    elements = []
    for element in single_elements:
        for shape in element["shapes"]:
            for order in element["orders"]:
                elements.append(FiniteElement(element["family"], shape, order))

    # Add the mixed elements
    elements += mixed_elements
    num_elements = len(elements)

    # Test all elements
    num_tests = 0
    msg = "Verifying evaluate_basis and evaluate_basis_derivatives for elements"
    info("\n" + msg + "\n" + len(msg)*"-")
    for i, ufl_element in enumerate(elements):
        num_tests += verify_element(num_elements, i + 1, ufl_element)

    # print results
    error = print_results(num_tests, ffc_fail, gcc_fail, run_fail, dif_cri, dif_acc, correct)

    if not error:
        # Remove temporary directory
        os.chdir(os.pardir)
        shutil.rmtree("tmp")
    pop_level()

    return error

if __name__ == "__main__":
    # sys.exit(main(sys.argv[1:]))
    sys.exit(main(INFO))
