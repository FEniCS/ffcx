#!/usr/bin/env python
"""
Tests of utilities for handling shapes, indice, strides etc.
"""

from six.moves import xrange
from uflacs.analysis.indexing import (shape_to_strides, multiindex_to_component,
                                        component_to_multiindex, indexing_to_component)

from operator import eq as equal

def test_shape_to_strides():
    assert () == shape_to_strides(())
    assert equal((1,),     shape_to_strides((3,)))
    assert equal((2,1),    shape_to_strides((3,2)))
    assert equal((4,1),    shape_to_strides((3,4)))
    assert equal((12,4,1), shape_to_strides((6,3,4)))

def test_multiindex_to_component_to_multiindex():
    sh = (2, 3, 5)
    strides = shape_to_strides(sh)
    for i in xrange(sh[2]):
        for j in xrange(sh[1]):
            for k in xrange(sh[0]):
                index = (k, j, i)
                c = multiindex_to_component(index, strides)
                index2 = component_to_multiindex(c, strides)
                assert index == index2

def test_indexing_to_component():
    assert equal(0, indexing_to_component(  (), (),   ()))
    assert equal(0, indexing_to_component((0,), (), (2,)))
    assert equal(1, indexing_to_component((1,), (), (2,)))
    assert equal(3, indexing_to_component((1,1), (), (2,2)))
    for i in xrange(5):
        for j in xrange(3):
            for k in xrange(2):
                assert equal(15*k+5*j+i, indexing_to_component((k,j,i), (), (2,3,5)))
    # TODO: Add free indices to the mix!
