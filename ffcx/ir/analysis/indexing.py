# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Algorithms for working with multiindices."""

import ufl
from ufl.classes import ComponentTensor, FixedIndex, Index, Indexed
from ufl.permutation import compute_indices
from ufl.utils.indexflattening import flatten_multiindex, shape_to_strides


def map_indexed_arg_components(indexed):
    """Build a map from flattened components to subexpression.

    Builds integer list mapping between flattened components
    of indexed expression and its underlying tensor-valued subexpression.
    """
    assert isinstance(indexed, Indexed)

    # AKA indexed = tensor[multiindex]
    tensor, multiindex = indexed.ufl_operands

    # AKA e1 = e2[multiindex]
    # (this renaming is historical, but kept for consistency with all the variables *1,*2 below)
    e2 = tensor
    e1 = indexed

    # Get tensor and index shape
    sh1 = e1.ufl_shape
    sh2 = e2.ufl_shape
    fi1 = e1.ufl_free_indices
    fi2 = e2.ufl_free_indices
    fid1 = e1.ufl_index_dimensions
    fid2 = e2.ufl_index_dimensions

    # Compute regular and total shape
    tsh1 = sh1 + fid1
    tsh2 = sh2 + fid2
    # r1 = len(tsh1)
    r2 = len(tsh2)
    # str1 = shape_to_strides(tsh1)
    str2 = shape_to_strides(tsh2)
    assert not sh1
    assert sh2  # Must have shape to be indexed in the first place
    assert ufl.product(tsh1) <= ufl.product(tsh2)

    # Build map from fi2/fid2 position (-offset nmui) to fi1/fid1 position
    ind2_to_ind1_map = [None] * len(fi2)
    for k, i in enumerate(fi2):
        ind2_to_ind1_map[k] = fi1.index(i)

    # Build map from fi1/fid1 position to mi position
    nmui = len(multiindex)
    multiindex_to_ind1_map = [None] * nmui
    for k, i in enumerate(multiindex):
        if isinstance(i, Index):
            multiindex_to_ind1_map[k] = fi1.index(i.count())

    # Build map from flattened e1 component to flattened e2 component
    perm1 = compute_indices(tsh1)
    ni1 = ufl.product(tsh1)

    # Situation: e1 = e2[mi]
    d1 = [None] * ni1
    p2 = [None] * r2
    assert len(sh2) == nmui
    for k, i in enumerate(multiindex):
        if isinstance(i, FixedIndex):
            p2[k] = int(i)
    for c1, p1 in enumerate(perm1):
        for k, i in enumerate(multiindex):
            if isinstance(i, Index):
                p2[k] = p1[multiindex_to_ind1_map[k]]
        for k, i in enumerate(ind2_to_ind1_map):
            p2[nmui + k] = p1[i]
        c2 = flatten_multiindex(p2, str2)
        d1[c1] = c2

    # Consistency checks
    assert all(isinstance(x, int) for x in d1)
    assert len(set(d1)) == len(d1)
    return d1


def map_component_tensor_arg_components(tensor):
    """Build a map from flattened components to subexpression.

    Builds integer list mapping between flattended components
    of tensor and its underlying indexed subexpression.
    """
    assert isinstance(tensor, ComponentTensor)

    # AKA tensor = as_tensor(indexed, multiindex)
    indexed, multiindex = tensor.ufl_operands

    e1 = indexed
    e2 = tensor  # e2 = as_tensor(e1, multiindex)
    mi = [i for i in multiindex if isinstance(i, Index)]

    # Get tensor and index shapes
    sh1 = e1.ufl_shape  # (sh)ape of e1
    sh2 = e2.ufl_shape  # (sh)ape of e2
    fi1 = e1.ufl_free_indices  # (f)ree (i)ndices of e1
    fi2 = e2.ufl_free_indices  # ...
    fid1 = e1.ufl_index_dimensions  # (f)ree (i)ndex (d)imensions of e1
    fid2 = e2.ufl_index_dimensions  # ...

    # Compute total shape (tsh) of e1 and e2
    tsh1 = sh1 + fid1
    tsh2 = sh2 + fid2
    r1 = len(tsh1)  # 'total rank' or e1
    r2 = len(tsh2)  # ...
    str1 = shape_to_strides(tsh1)
    assert not sh1
    assert sh2
    assert len(mi) == len(multiindex)
    assert ufl.product(tsh1) == ufl.product(tsh2)
    assert fi1

    assert all(i in fi1 for i in fi2)

    nmui = len(multiindex)
    assert nmui == len(sh2)

    # Build map from fi2/fid2 position (-offset nmui) to fi1/fid1 position
    p2_to_p1_map = [None] * r2
    for k, i in enumerate(fi2):
        p2_to_p1_map[k + nmui] = fi1.index(i)

    # Build map from fi1/fid1 position to mi position
    for k, i in enumerate(mi):
        p2_to_p1_map[k] = fi1.index(mi[k].count())

    # Build map from flattened e1 component to flattened e2 component
    perm2 = compute_indices(tsh2)
    ni2 = ufl.product(tsh2)

    # Situation: e2 = as_tensor(e1, mi)
    d2 = [None] * ni2
    p1 = [None] * r1
    for c2, p2 in enumerate(perm2):
        for k2, k1 in enumerate(p2_to_p1_map):
            p1[k1] = p2[k2]
        c1 = flatten_multiindex(p1, str1)
        d2[c2] = c1

    # Consistency checks
    assert all(isinstance(x, int) for x in d2)
    assert len(set(d2)) == len(d2)
    return d2
