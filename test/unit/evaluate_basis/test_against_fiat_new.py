# Copyright (C) 2010-2016 Kristian B. Oelgaard and Garth N. Wells
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

import os
import numpy
import time
import pytest

from ufl import FiniteElement, MixedElement

from ffc.log import push_level, pop_level, CRITICAL, INFO
from ffc.log import info, error, debug, info_blue, info_green
from ffc.mixedelement import MixedElement as FFCMixedElement

from ffc.fiatinterface import create_element, reference_cell

from instant.output import get_status_output

# Local imports
from cppcode import evaluate_basis_code_fiat
from test_common import compile_gcc_code, run_code, xcomb, get_element_name,\
    compile_element

tol = 1e-14
crit_tol = 1e-8

ffc_fail = []
gcc_fail = []
run_fail = []
dif_cri = []
dif_acc = []
correct = []

log_file = "fiat_errors.log"

# Some random points
random_points = {1: [(0.114,), (0.349,), (0.986,)],
                 2: [(0.114, 0.854), (0.349, 0.247), (0.986, 0.045)],
                 3: [(0.114, 0.854, 0.126), (0.349, 0.247, 0.457),
                     (0.986, 0.045, 0.127)]}


# Elements, supported by FFC and FIAT, and their supported shape and orders
single_elements = [{"family": "Lagrange",
                    "shapes": ["interval", "triangle", "tetrahedron"],
                    "orders": [1, 2, 3, 4]},
                   {"family": "Discontinuous Lagrange",
                    "shapes": ["interval", "triangle", "tetrahedron"],
                    "orders": [0, 1, 2, 3, 4]},
                   {"family": "Crouzeix-Raviart",
                    "shapes": ["triangle", "tetrahedron"],
                    "orders": [1]},
                   {"family": "Raviart-Thomas",
                    "shapes": ["triangle", "tetrahedron"],
                    "orders": [1, 2, 3]},
                   {"family": "Discontinuous Raviart-Thomas",
                    "shapes": ["triangle", "tetrahedron"],
                    "orders": [1, 2, 3]},
                   {"family": "Brezzi-Douglas-Marini",
                    "shapes": ["triangle", "tetrahedron"],
                    "orders": [1, 2, 3]},
                   {"family": "Brezzi-Douglas-Fortin-Marini",
                    "shapes": ["triangle"],
                    "orders": [2]},
                   {"family": "Nedelec 1st kind H(curl)",
                    "shapes": ["triangle", "tetrahedron"],
                    "orders": [1, 2, 3]},
                   {"family": "Nedelec 2nd kind H(curl)",
                    "shapes": ["triangle", "tetrahedron"],
                    "orders": [1, 2, 3]},
                   {"family": "Regge",
                    "shapes": ["triangle", "tetrahedron"],
                    "orders": [0, 1, 2, 3]}]

# import json
# with open('data.txt', 'w') as outfile:
#    json.dump(single_elements, outfile, indent=4)

# Create some mixed elements
dg0_tri = FiniteElement("DG", "triangle", 0)
dg1_tri = FiniteElement("DG", "triangle", 1)
cg1_tri = FiniteElement("CG", "triangle", 1)
cr1_tri = FiniteElement("CR", "triangle", 1)
rt1_tri = FiniteElement("RT", "triangle", 1)
drt2_tri = FiniteElement("DRT", "triangle", 2)
bdm1_tri = FiniteElement("BDM", "triangle", 1)
ned1_tri = FiniteElement("N1curl", "triangle", 1)
reg0_tri = FiniteElement("Regge", "triangle", 0)

dg0_tet = FiniteElement("DG", "tetrahedron", 0)
dg1_tet = FiniteElement("DG", "tetrahedron", 1)
cg1_tet = FiniteElement("CG", "tetrahedron", 1)
cr1_tet = FiniteElement("CR", "tetrahedron", 1)
rt1_tet = FiniteElement("RT", "tetrahedron", 1)
drt2_tet = FiniteElement("DRT", "tetrahedron", 2)
bdm1_tet = FiniteElement("BDM", "tetrahedron", 1)
ned1_tet = FiniteElement("N1curl", "tetrahedron", 1)
reg0_tet = FiniteElement("Regge", "tetrahedron", 0)


def matrix(points):
    return "{%s};" % ", ".join(["{%s}" % ", ".join([str(c) for c in p]) for p in points])


# Create combinations in pairs
mix_tri = [MixedElement(e) for e in xcomb([dg0_tri, dg1_tri, cg1_tri, cr1_tri,
                                           rt1_tri, drt2_tri, bdm1_tri,
                                           ned1_tri, reg0_tri], 2)]
mix_tet = [MixedElement(e) for e in xcomb([dg0_tet, dg1_tet, cg1_tet, cr1_tet,
                                           rt1_tet, drt2_tet, bdm1_tet,
                                           ned1_tet, reg0_tet], 2)]

mixed_elements = [MixedElement([dg0_tri]*4),
                  MixedElement([cg1_tri]*3),
                  MixedElement([bdm1_tri]*2),
                  MixedElement([dg1_tri, cg1_tri, cr1_tri, rt1_tri, bdm1_tri,
                                ned1_tri]),
                  MixedElement([MixedElement([rt1_tri, cr1_tri]), cg1_tri,
                                ned1_tri]),
                  MixedElement([ned1_tri, dg1_tri, MixedElement([rt1_tri,
                                                                 cr1_tri])]),
                  MixedElement([drt2_tri, cg1_tri]),
                  MixedElement([dg0_tet]*4), MixedElement([cg1_tet]*3),
                  MixedElement([bdm1_tet]*2),
                  MixedElement([dg1_tet, cg1_tet, cr1_tet, rt1_tet, bdm1_tet,
                                ned1_tet]),
                  MixedElement([MixedElement([rt1_tet, cr1_tet]), cg1_tet,
                                ned1_tet]),
                  MixedElement([ned1_tet, dg1_tet, MixedElement([rt1_tet,
                                                                 cr1_tet])]),
                  MixedElement([drt2_tet, cg1_tet]),
                  MixedElement([cg1_tet, cg1_tet, cg1_tet, reg0_tet])] + mix_tri + mix_tet


def to_fiat_tuple(comb, geo_dim):
    """Convert a list of combinations of derivatives to a fiat tuple of
    derivatives.  FIAT expects a list with the number of derivatives
    in each spatial direction.  E.g., in 2D: u_{xyy} --> [0, 1, 1] in
    FFC --> (1, 2) in FIAT.

    """
    new_comb = [0]*geo_dim
    if comb == []:
        return tuple(new_comb)

    for i in range(geo_dim):
        new_comb[i] = comb.count(i)

    return tuple(new_comb)


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


def get_data(ufl_element):
    "Get needed data to run tests."

    # Create fiat element.
    element = create_element(ufl_element)

    # The derivative order that we are interested in is the degree of
    # the element.
    if isinstance(element, FFCMixedElement):
        deriv_order = max([e.degree() for e in element.elements()])
    else:
        deriv_order = element.degree()

    # Get coordinates of the reference cell.
    cell = ufl_element.cell()
    ref_coords = reference_cell(cell.cellname()).get_vertices()

    # Get the locations of the fiat element dofs.
    elem_points = [list(L.pt_dict.keys())[0] for L in element.dual_basis()]

    # Add some random points.
    geo_dim = cell.geometric_dimension()
    points = elem_points + random_points[geo_dim]

    return (element, points, geo_dim, ref_coords, deriv_order)


def get_ffc_values(ufl_element):
    "Get the values from evaluate_basis and evaluate_basis_derivatives."

    # Get data from element and tabulate basis values
    element, points, geo_dim, ref_coords, deriv_order = get_data(ufl_element)

    # Get relevant element name
    element_name = get_element_name(ufl_element)

    # Create g++ code and compile
    num_coords = len(ref_coords)
    options = {"element": element_name,
               "dim": geo_dim,
               "num_points": len(points),
               "points": matrix(points),
               "cell_ref_coords": "double cell_ref_coords[{}][{}] = {}".format(num_coords, geo_dim, matrix(ref_coords)),
               "num_coords": num_coords}
    error = compile_gcc_code(ufl_element, evaluate_basis_code_fiat % options,
                             gcc_fail, log_file)
    print(error)
    assert not error

    # Loop over derivative order and compute values
    ffc_values = {}
    for n in range(deriv_order + 1):
        values = run_code(ufl_element, n, run_fail, log_file)
        assert values is not None
        ffc_values[n] = values

    return ffc_values


def get_fiat_values(ufl_element):
    """Create a FIAT element and use it to tabulate the values on the
    reference element. The return values is a dictionary with keys
    equal to the derivative order and values is a matrix where each
    row is the basis values at a point.  E.g., {0:[[1,0,0],[0,1,0],
    [0,0,1]]}.

    """

    # Get data and tabulate basis values.
    element, points, geo_dim, ref_coords, deriv_order = get_data(ufl_element)
    values = element.tabulate(deriv_order, points)
    return_values = {}
    value_shape = element.value_shape()

    # Rearrange values to match what we get from evaluate_basis*()
    for n in range(deriv_order + 1):
        combinations = compute_derivative_combinations(n, geo_dim)
        vals = []
        # If we're evaluating the basis functions, use all points, but
        # if we're evaluating the derivatives, just use the 3
        # arbitrary points to avoid the number of tests exploding with
        # spacedim**2.
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


def verify_values(ufl_element, ref_values, ffc_values, dif_cri, dif_acc, correct,
                  log_file):
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
            if error > crit_tol:
                m = "%s failed: error = %s (crit_tol: %s)" % (s, str(error), str(crit_tol))
                info_red(m)
                dif_cri.append(str(ufl_element))
                s = s + "\n" + m
                raise RuntimeError("%s failed: error = %s (crit_tol: %s)" % (s, str(error), str(crit_tol)))
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
    #if len(errors) > 1:
    #    log_error("\n".join(errors), log_file)

    return num_tests


@pytest.mark.parametrize("ufl_element", mixed_elements)
def test_element(ufl_element):
    "Test FFC elements against FIAT"

    info("\nVerifying element: {}".format(str(ufl_element)))

    # Compile element
    error = compile_element(ufl_element, ffc_fail, log_file)
    assert not error

    # Get FFC values
    t = time.time()
    ffc_values = get_ffc_values(ufl_element)
    if ffc_values is None:
        return
    debug("  time to compute FFC values: %f" % (time.time() - t))

    # Get FIAT values that are formatted in the same way as the values
    # from evaluate_basis and evaluate_basis_derivatives.
    fiat_values = get_fiat_values(ufl_element)

    # Compare FIAT and FFC values
    verify_values(ufl_element, fiat_values, ffc_values, dif_cri,
                  dif_acc, correct, log_file)
