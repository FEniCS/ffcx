import itertools
import numpy as np
import pytest

import ffc
import ffc_test_factory
import ffc_test_factory.factory


@pytest.fixture()
def build_element_list():
    element_data = ffc_test_factory.load()
    elements = [(ufl_e, i) for family, data in element_data.items() for i, name, ufl_e, index in data]
    element_ids = [str(ufl_e) for family, data in element_data.items() for i, name, ufl_e, index in data]

    return elements, element_ids


@pytest.fixture()
def point_data():
    # Some random points (1D, 3D, 3D)
    points = {1: [(0.114,), (0.349,), (0.986,)],
              2: [(0.114, 0.854), (0.349, 0.247), (0.986, 0.045)],
              #2: [(0.114, 0.854)],
              3: [(0.114, 0.854, 0.126), (0.349, 0.247, 0.457),
                  (0.986, 0.045, 0.127)]}

    return points


# @pytest.fixture(params=elements, ids=element_ids)
# @pytest.fixture(params=elements, ids=build_element_list)
@pytest.fixture(params=build_element_list()[0],
                ids=build_element_list()[1])
def element(request, point_data):

    # Extract UFC index and UFL element
    ufl_element, i = request.param

    # Create UFC element
    ufc_element = ffc_test_factory.factory.create_element(i)

    # Get geometric dim
    gdim = ufc_element.geometric_dimension()

    return ufl_element, ufc_element, point_data[gdim]


def test_evaluate_reference_basis_vs_fiat(element):
    """Tests ufc::finite_element::evaluate_reference_basis against data
    from FIAT

    """
    ufl_element, ufc_element, points = element

    # Get geometric and topological dimensions
    tdim = ufc_element.topological_dimension()

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
    # print("UFC __________   ", values_ufc.shape)
    assert np.allclose(values_ufc, values_fiat)

    # FIXME: This could be fragile because it depend on the order of
    # the sub-element in UFL being the same at the UFC order.

    # Test sub-elements recursively
    n = ufc_element.num_sub_elements()
    ufl_subelements = ufl_element.sub_elements()
    assert n == len(ufl_subelements)
    for i, ufl_e in enumerate(ufl_subelements):
        ufc_sub_element = ufc_element.create_sub_element(i)
        test_evaluate_reference_basis_vs_fiat((ufl_e, ufc_sub_element, points))


def test_evaluate_basis_vs_fiat(element):
    """Tests ufc::finite_element::evaluate_basis and
    ufc::finite_element::evaluate_basis_all against data from FIAT

    """
    ufl_element, ufc_element, points = element

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
        test_evaluate_basis_vs_fiat((ufl_e, ufc_sub_element, points))


@pytest.mark.parametrize("order", range(4))
def test_evaluate_reference_basis_deriv_vs_fiat(element, order):
    """Tests ufc::finite_element::evaluate_reference_basis_derivatives against data
    from FIAT

    """

    ufl_element, ufc_element, points = element

    # Get geometric and topological dimensions
    tdim = ufc_element.geometric_dimension()
    gdim = ufc_element.topological_dimension()

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
        test_evaluate_reference_basis_deriv_vs_fiat((ufl_e, ufc_sub_element, points), order)
