# Copyright (C) 2004-2009 Anders Logg
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
# Modified by Garth N. Wells, 2006
# Modified by Marie E. Rognes, 2007
# Modified by Kristian B. Oelgaard, 2009
#
# First added:  2004-11-03
# Last changed: 2009-12-21

# Python modules.
import functools
import numpy

# FFC modules.
from ffc.utils import listcopy
from ffc.log import error

def build_indices(dims):
    "Create a list of all index combinations."
    if not dims: return [[]]
    ranges = listcopy(dims)
    return functools.reduce(outer_join, ranges, [[]])

def outer_join(a, b):
    """Let a be a list of lists and b a list. We append each element
    of b to each list in a and return the resulting list of lists."""
    outer = []
    for i in range(len(a)):
        for j in range(len(b)):
            outer += [a[i] + [b[j]]]
    return outer

def create_multiindex(indices):
    "Create multiindex for given list of indices."

    # Check that we got all indices correctly
    indices = sorted(indices)
    for (i, index) in enumerate(indices):
        if not i == index.index_id:
            error("Unable to extract all indices.")

    # Get dimensions
    dims = [range(len(index.index_range)) for index in indices]

    return MultiIndex(dims)

class MultiIndex:
    """
    A MultiIndex represents a list of indices and holds the following
    data:

        rank    - rank of multiindex
        dims    - a list of dimensions
        indices - a list of all possible multiindex values
    """

    def __init__(self, dims):
        "Create multiindex from given list of ranges"
        self.rank = len(dims)
        self.dims = [len(dim) for dim in dims]
        self.indices = build_indices(dims)
        return

    def __str__(self):
        "Return informal string representation (pretty-print)."
        return "rank = %d dims = %s indices = %s" % (self.rank, str(self.dims), str(self.indices))
