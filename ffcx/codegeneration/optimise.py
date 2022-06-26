# Copyright (C) 2022 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections

from ffcx.codegeneration.indices import MultiIndex


def fuse_loops(lang, definitions):
    """
    Merge a sequence of loops with the same iteration space into a single loop.

    Loop fusion improves data locality, cache reuse and decreases the loop control overhead.
    NOTE: Loop fusion might increase the pressure on register allocation.
    Ideally, we should define a cost function to determine how many loops should fuse at a time.
    """

    bodies = collections.defaultdict(list)
    indices = collections.defaultdict(MultiIndex)

    pre_loop = []

    for access, definition in definitions.items():
        for d in definition:
            if isinstance(d, lang.NestedForRange):
                index = d.multi_indices[0]
                if index.dim == 1:
                    hash_ = hash(index)
                    bodies[hash_] += [d.body()]
                    indices[hash_] = index
                else:
                    pre_loop += [d]
            else:
                pre_loop += [d]

    fused = []
    for key in indices.keys():
        body = bodies[key]
        index = indices[key]
        fused += [lang.NestedForRange([index], body)]


    code = []
    code += pre_loop
    code += fused
    return code
