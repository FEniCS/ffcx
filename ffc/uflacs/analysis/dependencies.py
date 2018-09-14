# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tools for analysing dependencies within expression graphs."""

import numpy


def compute_dependencies(e2i, V, ignore_terminal_modifiers=True):

    dependencies = []
    for v in V:
        if v._ufl_is_terminal_ or (ignore_terminal_modifiers
                                   and v._ufl_is_terminal_modifier_):
            dependencies.append(())
        else:
            dependencies.append([e2i[o] for o in v.ufl_operands])

    return dependencies


def invert_dependencies(dependencies):
    n = len(dependencies)
    invdeps = [()] * n
    for i in range(n):
        for d in dependencies[i]:
            invdeps[d] += (i, )
    return invdeps


def mark_active(dependencies, targets):
    """Return an array marking the recursive dependencies of targets.

    Input:
    - dependencies - List of ints, a mapping from a symbol to the symbols of its dependencies.
    - targets      - Sequence of symbols to mark the dependencies of.

    Output:
    - active   - Truth value for each symbol.
    - num_used - Number of true values in active array.
    """
    n = len(dependencies)

    # Initial state where nothing is marked as used
    active = numpy.zeros(n, dtype=numpy.int8)
    num_used = 0

    # Seed with initially used symbols
    assert isinstance(targets, list)
    active[targets] = 1

    # Mark dependencies by looping backwards through symbols array
    for s in range(n - 1, -1, -1):
        if active[s]:
            num_used += 1
            active[list(dependencies[s])] = 1

    # Return array marking which symbols are used and the number of positives
    return active, num_used


def mark_image(inverse_dependencies, sources):
    """Return an array marking the set of symbols dependent on the sources.

    Input:
    - dependencies - List of ints, a mapping from a symbol to the symbols of its dependencies.
    - sources      - Sequence of symbols to mark the dependants of.

    Output:
    - image    - Truth value for each symbol.
    - num_used - Number of true values in active array.
    """
    n = len(inverse_dependencies)

    # Initial state where nothing is marked as used
    image = numpy.zeros(n, dtype=numpy.int8)
    num_used = 0

    # Seed with initially used symbols
    assert isinstance(sources, list)
    image[sources] = 1

    # Mark dependencies by looping forwards through symbols array
    for s in range(n):
        if image[s]:
            num_used += 1
            image[list(inverse_dependencies[s])] = 1

    # Return array marking which symbols are used and the number of positives
    return image, num_used
