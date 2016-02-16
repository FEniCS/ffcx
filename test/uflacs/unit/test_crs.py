#!/usr/bin/env python
"""
Tests of CRS data structure.
"""

from six.moves import xrange
from uflacs.datastructures.crs import CRS

def test_crs_can_have_zero_element_rows():
    rcap, ecap = 3, 1
    A = CRS(rcap, ecap, int)
    for i in range(rcap):
        row = []
        A.push_row(row)
    assert A.num_elements == 0
    for i in range(rcap):
        row = []
        assert list(A[i]) == row

def test_crs_can_have_one_element_rows():
    rcap, ecap = 3, 3
    A = CRS(rcap, ecap, int)
    for i in range(rcap):
        row = [i]
        A.push_row(row)
    assert A.num_elements == rcap
    for i in range(rcap):
        row = [i]
        assert list(A[i]) == row

def test_crs_can_have_n_element_rows():
    rcap, ecap = 5, 25
    A = CRS(rcap, ecap, int)
    k = 0
    for i in range(rcap):
        row = [i+2, i+1] + [i]*i
        k += len(row)
        A.push_row(row)
    assert A.num_elements == k
    for i in range(rcap):
        row = [i+2, i+1] + [i]*i
        assert list(A[i]) == row
