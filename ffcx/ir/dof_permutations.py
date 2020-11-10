# Copyright (C) 2020 Matthew W. Scroggs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import warnings
import math
import ufl
import numpy as np


def block_size_perm(perm, block_size):
    new_perm = np.zeros((perm.shape[0] * block_size, perm.shape[1] * block_size))
    for i, row in enumerate(perm):
        for j, entry in enumerate(row):
            new_perm[i * block_size: (i+1) * block_size, j * block_size: (j+1) * block_size] = entry
    return new_perm


def blocked_perm(perms):
    new_perm = np.zeros((sum(i.shape[0] for i in perms), sum(i.shape[1] for i in perms)))
    row_start = 0
    col_start = 0
    for i in perms:
        new_perm[row_start: row_start + i.shape[0], col_start: col_start + i.shape[1]] = i
        row_start += i.shape[0]
        col_start += i.shape[1]
    return new_perm


def base_permutations(ufl_element, elements):
    """Returns the base permutations."""
    if len(elements) == 1:
        return elements[0].base_permutations

    assert ufl_element.num_sub_elements() != 0

    # If the element has sub elements, combine their permutations
    perms = None

    if isinstance(ufl_element, ufl.VectorElement) or isinstance(ufl_element, ufl.TensorElement):
        block_size = ufl_element.num_sub_elements()
        return [
            block_size_perm(perm, block_size)
            for perm in base_permutations(ufl_element.sub_elements()[0], elements[:1])
        ]

    for e, e2 in zip(ufl_element.sub_elements(), elements):
        bp = base_permutations(e, [e2])
        if perms is None:
            perms = [[] for i in bp]
        for i, b in enumerate(bp):
            perms[i] += [a + len(perms[i]) for a in b]
    return [blocked_perm(p) for p in perms]
