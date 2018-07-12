# -*- coding: utf-8 -*-
# Copyright (C) 2017 Garth N. Wells
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

import numpy as np
import pytest

from ufl import FiniteElement, VectorElement, MixedElement
from ufl import interval, triangle, tetrahedron
import ffc
import ffc_factory


@pytest.fixture(scope="module")
def point_data():
    """Points are which evaluate functions are tested"""
    p = {1: [(0.114,), (0.349,), (0.986,)],
         2: [(0.114, 0.854), (0.349, 0.247), (0.986, 0.045)],
         3: [(0.114, 0.854, 0.126), (0.349, 0.247, 0.457),
             (0.986, 0.045, 0.127)]}
    return p


@pytest.fixture(scope="module")
def build_ufl_element_list():
    """Build collection of UFL elements"""

    elements = []

    # Lagrange elements
    for cell in (interval, triangle, tetrahedron):
        for p in range(1, 2):
            elements.append(FiniteElement("Lagrange", cell, p))
            elements.append(VectorElement("Lagrange", cell, p))
            elements.append(FiniteElement("Discontinuous Lagrange", cell, p-1))

    # Vector elements
    for cell in (triangle, tetrahedron):
        for p in range(1, 2):
            elements.append(FiniteElement("RT", cell, p))
            elements.append(FiniteElement("BDM", cell, p))
            elements.append(FiniteElement("N1curl", cell, p))
            elements.append(FiniteElement("N2curl", cell, p))
            elements.append(FiniteElement("Discontinuous Raviart-Thomas", cell, p))

    # Mixed elements
    for cell in (interval, triangle, tetrahedron):
        for p in range(1, 3):
            e0 = FiniteElement("Lagrange", cell, p+1)
            e1 = FiniteElement("Lagrange", cell, p)
            e2 = VectorElement("Lagrange", cell, p+1)

            elements.append(MixedElement([e0, e0]))
            elements.append(MixedElement([e0, e1]))
            elements.append(MixedElement([e1, e0]))
            elements.append(MixedElement([e2, e1]))
            elements.append(MixedElement([MixedElement([e1, e1]), e0]))

    for cell in (triangle, tetrahedron):
        for p in range(1, 2):
            e0 = FiniteElement("Lagrange", cell, p+1)
            e1 = FiniteElement("Lagrange", cell, p)
            e2 = VectorElement("Lagrange", cell, p+1)
            e3 = FiniteElement("BDM", cell, p)
            e4 = FiniteElement("RT", cell, p+1)
            e5 = FiniteElement("N1curl", cell, p)

            elements.append(MixedElement([e1, e2]))
            elements.append(MixedElement([e3, MixedElement([e4, e5])]))

    # Misc elements
    for cell in (triangle,):
        for p in range(1, 2):
            elements.append(FiniteElement("HHJ", cell, p))

    for cell in (triangle, tetrahedron):
        elements.append(FiniteElement("CR", cell, 1))
        for p in range(1, 2):
            elements.append(FiniteElement("Regge", cell, p))

    return elements


def element_id(e):
    """Return element string signature as ID"""
    return str(e)


@pytest.fixture(params=build_ufl_element_list(), ids=element_id, scope="module")
def element_pair(request):
    """Given a UFL element, returns UFL element and a JIT-compiled and
    wrapped UFC element.

    """

    ufl_element = request.param
    ufc_element, ufc_dofmap = ffc.jit(ufl_element, parameters=None)
    ufc_element = ffc_factory.make_ufc_finite_element(ufc_element)
    return ufl_element, ufc_element


def test_evaluate_reference_basis_vs_fiat(element_pair, point_data):
    """Tests ufc::finite_element::evaluate_reference_basis against data
    from FIAT.

    """

    ufl_element, ufc_element = element_pair

    # Get geometric and topological dimensions
    tdim = ufc_element.topological_dimension()

    points = point_data[tdim]

    # Create FIAT element via FFC
    fiat_element = ffc.fiatinterface.create_element(ufl_element)

    # Tabulate basis at requsted points via UFC
    values_ufc = ufc_element.evaluate_reference_basis(points)

    # Tabulate basis at requsted points via FIAT
    values_fiat = fiat_element.tabulate(0, points)

    # Extract and reshape FIAT output
    key = (0,)*tdim
    values_fiat = np.array(values_fiat[key])

    # Shape checks
    fiat_value_size = 1
    for i in range(ufc_element.reference_value_rank()):
        fiat_value_size *= values_fiat.shape[i + 1]
        assert values_fiat.shape[i + 1] == ufc_element.reference_value_dimension(i)

    # Reshape FIAT output to compare with UFC output
    values_fiat = values_fiat.reshape((values_fiat.shape[0],
                                       ufc_element.reference_value_size(),
                                       values_fiat.shape[-1]))

    # Transpose
    values_fiat = values_fiat.transpose((2, 0, 1))

    # Check values
    assert np.allclose(values_ufc, values_fiat)

    # FIXME: This could be fragile because it depend on the order of
    # the sub-element in UFL being the same at the UFC order.

    # Test sub-elements recursively
    n = ufc_element.num_sub_elements()
    ufl_subelements = ufl_element.sub_elements()
    assert n == len(ufl_subelements)
    for i, ufl_e in enumerate(ufl_subelements):
        ufc_sub_element = ufc_element.create_sub_element(i)
        test_evaluate_reference_basis_vs_fiat((ufl_e, ufc_sub_element), point_data)


def test_evaluate_basis_vs_fiat(element_pair, point_data):
    """Tests ufc::finite_element::evaluate_basis and
    ufc::finite_element::evaluate_basis_all against data from FIAT.

    """

    ufl_element, ufc_element = element_pair

    # Get geometric and topological dimensions
    gdim = ufc_element.geometric_dimension()

    points = point_data[gdim]

    # Get geometric and topological dimensions
    tdim = ufc_element.topological_dimension()

    # Create FIAT element via FFC
    fiat_element = ffc.fiatinterface.create_element(ufl_element)

    # Tabulate basis at requsted points via FIAT
    values_fiat = fiat_element.tabulate(0, points)

    # Extract and reshape FIAT output
    key = (0,)*tdim
    values_fiat = np.array(values_fiat[key])

    # Shape checks
    fiat_value_size = 1
    for i in range(ufc_element.reference_value_rank()):
        fiat_value_size *= values_fiat.shape[i + 1]
        assert values_fiat.shape[i + 1] == ufc_element.reference_value_dimension(i)

    # Reshape FIAT output to compare with UFC output
    values_fiat = values_fiat.reshape((values_fiat.shape[0],
                                       ufc_element.reference_value_size(),
                                       values_fiat.shape[-1]))

    # Transpose
    values_fiat = values_fiat.transpose((2, 0, 1))

    # Get reference cell
    cell = ufl_element.cell()
    ref_coords = ffc.fiatinterface.reference_cell(cell.cellname()).get_vertices()

    # Iterate over each point and test
    for i, p in enumerate(points):
        values_ufc = ufc_element.evaluate_basis_all(p, ref_coords, 1)
        assert np.allclose(values_ufc, values_fiat[i])

        for d in range(ufc_element.space_dimension()):
            values_ufc = ufc_element.evaluate_basis(d, p, ref_coords, 1)
            assert np.allclose(values_ufc, values_fiat[i][d])

    # Test sub-elements recursively
    n = ufc_element.num_sub_elements()
    ufl_subelements = ufl_element.sub_elements()
    assert n == len(ufl_subelements)
    for i, ufl_e in enumerate(ufl_subelements):
        ufc_sub_element = ufc_element.create_sub_element(i)
        test_evaluate_basis_vs_fiat((ufl_e, ufc_sub_element), point_data)


@pytest.mark.parametrize("order", range(4))
def test_evaluate_reference_basis_deriv_vs_fiat(order, element_pair,
                                                point_data):
    """Tests ufc::finite_element::evaluate_reference_basis_derivatives
    against data from FIAT.

    """

    ufl_element, ufc_element = element_pair

    # Get geometric and topological dimensions
    tdim = ufc_element.geometric_dimension()
    gdim = ufc_element.topological_dimension()

    points = point_data[tdim]

    # Create FIAT element via FFC
    fiat_element = ffc.fiatinterface.create_element(ufl_element)

    # Tabulate basis at requsted points via UFC
    values_ufc = ufc_element.evaluate_reference_basis_derivatives(order, points)

    # Tabulate basis at requsted points via FIAT
    values_fiat = fiat_element.tabulate(order, points)

    # Compute 'FFC' style derivative indicators
    deriv_order = order
    geo_dim = tdim
    num_derivatives = geo_dim**deriv_order
    combinations = [[0] * deriv_order for n in range(num_derivatives)]
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

    def to_fiat_tuple(comb, gdim):
        """Convert a list of combinations of derivatives to a fiat tuple of
        derivatives.  FIAT expects a list with the number of
        derivatives in each spatial direction.  E.g., in 2D: u_{xyy}
        --> [0, 1, 1] in FFC --> (1, 2) in FIAT.

        """
        new_comb = [0] * gdim
        if comb == []:
            return tuple(new_comb)
        for i in range(gdim):
            new_comb[i] = comb.count(i)
        return tuple(new_comb)

    derivs = [to_fiat_tuple(c, tdim) for c in combinations]
    assert len(derivs) == tdim**order

    # Iterate over derivatives
    for i, d in enumerate(derivs):
        values_fiat_slice = np.array(values_fiat[d])

        # Shape checks
        fiat_value_size = 1
        for j in range(ufc_element.reference_value_rank()):
            fiat_value_size *= values_fiat_slice.shape[j + 1]
            assert values_fiat_slice.shape[j + 1] == ufc_element.reference_value_dimension(j)

        # Reshape FIAT output to compare with UFC output
        values_fiat_slice = values_fiat_slice.reshape((values_fiat_slice.shape[0],
                                                       ufc_element.reference_value_size(),
                                                       values_fiat_slice.shape[-1]))

        # Transpose
        values_fiat_slice = values_fiat_slice.transpose((2, 0, 1))

        # Compare
        assert np.allclose(values_ufc[:,:,i,:], values_fiat_slice)

    # Test sub-elements recursively
    n = ufc_element.num_sub_elements()
    ufl_subelements = ufl_element.sub_elements()
    assert n == len(ufl_subelements)
    for i, ufl_e in enumerate(ufl_subelements):
        ufc_sub_element = ufc_element.create_sub_element(i)
        test_evaluate_reference_basis_deriv_vs_fiat(order, (ufl_e, ufc_sub_element), point_data)
