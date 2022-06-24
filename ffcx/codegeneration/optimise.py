# Copyright (C) 2022 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
import collections


def fuse_loops(lang, definitions):
    """
    Merge a sequence of loops with the same iteration space into a single loop.

    Loop fusion improves data locality, cache reuse and decreases the loop control overhead.
    NOTE: Loop fusion might increase the pressure on register allocation.
    Ideally, we should define a cost function to determine how many loops should fuse at a time.
    """

    loops = collections.defaultdict(list)
    pre_loop = []
    for access, definition in definitions.items():
        
        for d in definition:
            if isinstance(d, lang.NestedForRange):
                index = d.multi_index
                loops[index].append(d.body())
                print(len(loops[index]))
            else:
                pre_loop += [d]
    fused = []


    for multi_index, body in loops.items():
        fused += [lang.NestedForRange(multi_index, body)]

    code = []
    code += pre_loop
    code += fused
    return code
