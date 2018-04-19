# -*- coding: utf-8 -*-
# Copyright (C) 2005-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffc import FFCError

logger = logging.getLogger(__name__)


def all_equal(sequence):
    """Check that all items in list are equal."""
    return sequence[:-1] == sequence[1:]


def pick_first(sequence):
    """Check that all values are equal and return the value."""
    if not all_equal(sequence):
        raise FFCError("Values differ: {}".format(sequence))
    return sequence[0]


def listcopy(sequence):
    """Create a copy of the list, calling the copy constructor on each
    object in the list (problems when using copy.deepcopy).

    """
    if not sequence:
        return []
    else:
        return [object.__class__(object) for object in sequence]


def compute_permutations(k, n, skip=[]):
    """Compute all permutations of k elements from (0, n) in rising order.
    Any elements that are contained in the list skip are not included."""
    if k == 1:
        return [(i, ) for i in range(n) if i not in skip]
    pp = compute_permutations(k - 1, n, skip)
    permutations = []
    for i in range(n):
        if i in skip:
            continue
        for p in pp:
            if i < p[0]:
                permutations += [(i, ) + p]

    return permutations


def insert_nested_dict(root, keys, value):
    """Set root[keys[0]][...][keys[-1]] = value, creating subdicts on the
    way if missing.

    """
    for k in keys[:-1]:
        d = root.get(k)
        if d is None:
            d = {}
            root[k] = d
        root = d
    root[keys[-1]] = value
