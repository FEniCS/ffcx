# -*- coding: utf-8 -*-
"""
Tests of table manipulation utilities.
"""

from ufl import triangle
from ffc.uflacs.elementtables import equal_tables, strip_table_zeros, build_unique_tables, get_ffc_table_values

import numpy as np
default_tolerance = 1e-14


def old_strip_table_zeros(table, eps):
    "Quick fix translating new version to old behaviour to keep tests without modification."
    dofrange, dofmap, stripped_table = strip_table_zeros(table, False, eps)
    begin, end = dofrange
    return begin, end, stripped_table


def test_equal_tables():
    a = np.zeros((2, 3))
    b = np.zeros((2, 3))
    assert equal_tables(a, b, default_tolerance)

    a = np.ones((2, 3))
    b = np.zeros((2, 3))
    assert not equal_tables(a, b, default_tolerance)

    a = np.ones((2, 3))
    b = np.ones((3, 2))
    assert not equal_tables(a, b, default_tolerance)

    a = np.ones((2, 3))*1.1
    b = np.ones((2, 3))
    assert not equal_tables(a, b, default_tolerance)

    # Checking for all ones within given tolerance:
    eps = 1e-4
    a = np.ones((2, 3))*(1.0+eps)
    assert equal_tables(a, np.ones(a.shape), 10*eps)


def test_strip_table_zeros():
    # Can strip entire table:
    a = np.zeros((2, 3))
    e = np.zeros((2, 0))
    begin, end, b = old_strip_table_zeros(a, default_tolerance)
    assert begin == end # This is a way to check for all-zero table
    assert begin == 0
    assert end == 0
    assert equal_tables(b, e, default_tolerance)

    # Can keep entire nonzero table:
    a = np.ones((2, 3))
    e = np.ones((2, 3))
    begin, end, b = old_strip_table_zeros(a, default_tolerance)
    assert begin != end
    assert begin == 0
    assert end == a.shape[-1]
    assert equal_tables(b, e, default_tolerance)

    # Can strip one left side column:
    a = np.ones((2, 3))
    a[:, 0] = 0.0
    e = np.ones((2, 2))
    begin, end, b = old_strip_table_zeros(a, default_tolerance)
    assert begin != end
    assert begin == 1
    assert end == a.shape[-1]
    assert equal_tables(b, e, default_tolerance)

    # Can strip one right side column:
    a = np.ones((2, 3))
    a[:, 2] = 0.0
    e = np.ones((2, 2))
    begin, end, b = old_strip_table_zeros(a, default_tolerance)
    assert begin != end
    assert begin == 0
    assert end == a.shape[-1]-1
    assert equal_tables(b, e, default_tolerance)

    # Can strip two columns on each side:
    a = np.ones((2, 5))
    a[:, 0] = 0.0
    a[:, 1] = 0.0
    a[:, 3] = 0.0
    a[:, 4] = 0.0
    e = np.ones((2, 1))
    begin, end, b = old_strip_table_zeros(a, default_tolerance)
    assert begin != end
    assert begin == 2
    assert end == a.shape[-1]-2
    assert equal_tables(b, e, default_tolerance)

    # Can strip two columns on each side of rank 1 table:
    a = np.ones((5,))
    a[..., 0] = 0.0
    a[..., 1] = 0.0
    a[..., 3] = 0.0
    a[..., 4] = 0.0
    e = np.ones((1,))
    begin, end, b = old_strip_table_zeros(a, default_tolerance)
    assert begin != end
    assert begin == 2
    assert end == a.shape[-1]-2
    assert equal_tables(b, e, default_tolerance)

    # Can strip two columns on each side of rank 3 table:
    a = np.ones((3, 2, 5))
    a[..., 0] = 0.0
    a[..., 1] = 0.0
    a[..., 3] = 0.0
    a[..., 4] = 0.0
    e = np.ones((3, 2, 1))
    begin, end, b = old_strip_table_zeros(a, default_tolerance)
    assert begin != end
    assert begin == 2
    assert end == a.shape[-1]-2
    assert equal_tables(b, e, default_tolerance)

    # Can strip columns on each side and compress the one in the middle:
    a = np.ones((2, 5))
    a[:, 0] = 0.0
    a[:, 1] = 1.0
    a[:, 2] = 0.0
    a[:, 3] = 3.0
    a[:, 4] = 0.0
    e = np.ones((2, 2))
    e[:, 0] = 1.0
    e[:, 1] = 3.0
    dofrange, dofmap, b = strip_table_zeros(a, True, default_tolerance)
    begin, end = dofrange
    assert begin != end
    assert begin == 1
    assert end == a.shape[-1] - 1
    assert list(dofmap) == [1, 3]
    assert equal_tables(b, e, default_tolerance)


def test_unique_tables_some_equal():
    tables = [
        np.zeros((2,)),
        np.ones((2,)),
        np.zeros((3,)),
        np.ones((2,)),
        np.ones((2,))*2,
        np.ones((2,))*2,
        ]
    unique, mapping = build_unique_tables(tables, default_tolerance)
    expected_unique = [
        np.zeros((2,)),
        np.ones((2,)),
        np.zeros((3,)),
        # np.ones((2,)),
        np.ones((2,))*2,
        # np.ones((2,))*2,
        ]
    expected_mapping = dict((i, v) for i, v in enumerate([0, 1, 2, 1, 3, 3]))
    assert mapping == expected_mapping
    assert len(set(mapping.values())) == len(unique)
    for i, t in enumerate(tables):
        assert equal_tables(t, unique[mapping[i]], default_tolerance)


def test_unique_tables_all_equal():
    tables = [np.ones((3, 5))*2.0]*6
    unique, mapping = build_unique_tables(tables, default_tolerance)
    expected_unique = [tables[0]]
    expected_mapping = dict((i, v) for i, v in enumerate([0]*6))
    assert mapping == expected_mapping
    assert len(set(mapping.values())) == len(unique)
    for i, t in enumerate(tables):
        assert equal_tables(t, unique[mapping[i]], default_tolerance)

def test_unique_tables_all_different():
    tables = [
        np.ones((2,)),
        np.ones((2, 3)),
        np.ones((2, 3, 4)),
        np.ones((2, 3, 4, 5)),
        ]
    unique, mapping = build_unique_tables(tables, default_tolerance)
    expected_unique = tables
    expected_mapping = dict((i, i) for i in range(len(tables)))
    assert mapping == expected_mapping
    assert len(set(mapping.values())) == len(unique)
    for i, t in enumerate(tables):
        assert equal_tables(t, unique[mapping[i]], default_tolerance)

def test_unique_tables_string_keys():
    tables = {
        'a': np.zeros((2,)),
        'b': np.ones((2,)),
        'c': np.zeros((3,)),
        'd': np.ones((2,)),
        'e': np.ones((2,))*2,
        'f': np.ones((2,))*2,
        }
    unique, mapping = build_unique_tables(tables, default_tolerance)
    expected_unique = [
        np.zeros((2,)),
        np.ones((2,)),
        np.zeros((3,)),
        # np.ones((2,)),
        np.ones((2,))*2,
        # np.ones((2,))*2,
        ]
    expected_mapping = { 'a': 0, 'b': 1, 'c': 2, 'd': 1, 'e': 3, 'f': 3 }
    assert mapping == expected_mapping
    assert len(set(mapping.values())) == len(unique)
    for i, t in tables.items():
        assert equal_tables(t, unique[mapping[i]], default_tolerance)


def xtest_get_ffc_table_values_scalar_cell():
    cell = triangle
    integral_type = "cell"
    entitytype = "cell"

    class MockElement:
        def value_shape(self): return ()

    element = MockElement()
    component = ()

    avg = None
    entity = 0
    for num_points in (1, 3):
        for num_dofs in (1, 5):
            arr = np.ones((num_dofs, num_points))
            for derivatives in [(0, 0), (1, 1)]:
                # Mocking table built by ffc
                ffc_tables = {
                    num_points: {
                        element: {
                            avg: { # avg
                                entity: { # entityid
                                    derivatives: arr
                                }
                            }
                        }
                    }
                }
                table = get_ffc_table_values(ffc_tables,
                    cell, integral_type,
                    element, avg,
                    entitytype, derivatives, component,
                    default_tolerance)
                assert equal_tables(table[0, ...], np.transpose(arr), default_tolerance)


def xtest_get_ffc_table_values_vector_facet():
    cell = triangle
    integral_type = "exterior_facet"
    entitytype = "facet"
    num_entities = 3

    class MockElement:
        def value_shape(self): return (2,)

    element = MockElement()
    num_components = 2

    avg = None
    for num_points in (1, 5):
        for num_dofs in (4, 7):
            # Make ones array of the right shape (all dimensions differ to detect algorithm bugs better)
            arr1 = np.ones((num_dofs, num_components, num_points))
            arrays = []
            for i in range(num_entities):
                arr = (i+1.0)*arr1 # Make first digit the facet number (1,2,3)
                for j in range(num_components):
                    arr[:, j,:] += 0.1*(j+1.0) # Make first decimal the component number (1,2)
                for j in range(num_dofs):
                    arr[j,:,:] += 0.01*(j+1.0) # Make second decimal the dof number
                for j in range(num_points):
                    arr[:,:, j] += 0.001*(j+1.0) # Make third decimal the point number
                arrays.append(arr)

            for derivatives in [(), (0,)]:
                # Mocking table built by ffc
                ffc_tables = {
                    num_points: {
                        element: {
                            avg: {
                                entity: { derivatives: arrays[entity] }
                                for entity in range(num_entities)
                            }
                        }
                    }
                }
                # Tables use flattened component, so we can loop over them as integers:
                for component in range(num_components):
                    table = get_ffc_table_values(ffc_tables,
                        cell, integral_type,
                        element, avg,
                        entitytype, derivatives, component,
                        default_tolerance)
                    for i in range(num_entities):
                        # print table[i,...]
                        assert equal_tables(table[i, ...], np.transpose(arrays[i][:, component,:]), default_tolerance)
