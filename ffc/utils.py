# -*- coding: utf-8 -*-

# Copyright (C) 2005-2017 Anders Logg
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
# Modified by Kristian B. Oelgaard, 2009
# Modified by Martin Sandve Alnæs 2014

# Python modules.
import operator
import functools
import itertools

# FFC modules.
from .log import error

from ufl.utils.sequences import product


def all_equal(sequence):
    "Check that all items in list are equal."
    return sequence[:-1] == sequence[1:]


def pick_first(sequence):
    "Check that all values are equal and return the value."
    if not all_equal(sequence):
        error("Values differ: " + str(sequence))
    return sequence[0]


def listcopy(sequence):
    """Create a copy of the list, calling the copy constructor on each
    object in the list (problems when using copy.deepcopy)."""
    if not sequence:
        return []
    else:
        return [object.__class__(object) for object in sequence]


def compute_permutations(k, n, skip=[]):
    """Compute all permutations of k elements from (0, n) in rising order.
    Any elements that are contained in the list skip are not included."""
    if k == 1:
        return [(i,) for i in range(n) if i not in skip]
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
    "Set root[keys[0]][...][keys[-1]] = value, creating subdicts on the way if missing."
    for k in keys[:-1]:
        d = root.get(k)
        if d is None:
            d = {}
            root[k] = d
        root = d
    root[keys[-1]] = value
