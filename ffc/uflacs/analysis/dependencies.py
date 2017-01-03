# -*- coding: utf-8 -*-
# Copyright (C) 2011-2016 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""Tools for analysing dependencies within expression graphs."""

import numpy

from six.moves import xrange as range

from ffc.uflacs.analysis.crsarray import CRSArray, sufficient_int


def compute_dependencies(e2i, V, ignore_terminal_modifiers=True):
    # Use numpy int type sufficient to hold num_rows
    num_rows = len(V)
    itype = sufficient_int(num_rows)

    # Preallocate CRSArray matrix of sufficient capacity
    num_nonzeros = sum(len(v.ufl_operands) for v in V)
    dependencies = CRSArray(num_rows, num_nonzeros, itype)
    for v in V:
        if v._ufl_is_terminal_ or (ignore_terminal_modifiers and v._ufl_is_terminal_modifier_):
            dependencies.push_row(())
        else:
            dependencies.push_row([e2i[o] for o in v.ufl_operands])

    return dependencies


def mark_active(dependencies, targets):
    """Return an array marking the recursive dependencies of targets.

    Input:
    - dependencies - CRSArray of ints, a mapping from a symbol to the symbols of its dependencies.
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
    active[targets] = 1

    # Mark dependencies by looping backwards through symbols array
    for s in range(n - 1, -1, -1):
        if active[s]:
            num_used += 1
            active[dependencies[s]] = 1

    # Return array marking which symbols are used and the number of positives
    return active, num_used


def mark_image(inverse_dependencies, sources):
    """Return an array marking the set of symbols dependent on the sources.

    Input:
    - dependencies - CRSArray of ints, a mapping from a symbol to the symbols of its dependencies.
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
    image[sources] = 1

    # Mark dependencies by looping forwards through symbols array
    for s in range(n):
        if image[s]:
            num_used += 1
            image[inverse_dependencies[s]] = 1

    # Return array marking which symbols are used and the number of positives
    return image, num_used
