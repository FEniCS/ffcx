# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
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

"""Tools for computing various shapes of ufl expressions.

The total shape is the regular shape tuple plus the index shape tuple.
The index shape tuple is the tuple of index dimensions of the free indices
of the expression, sorted by the count of the free indices.

The total shape of a tensor valued expression ``A`` and
``A[*indices(len(A.ufl_shape))]`` is therefore the same.
"""


def compute_index_shape(v):
    """Compute the 'index shape' of v."""
    return v.ufl_index_dimensions


def compute_all_shapes(v):
    """Compute the tensor-, index-, and total shape of an expr.

    Returns (shape, size, index_shape, index_size, total_shape, total_size).
    """
    shape = v.ufl_shape
    index_shape = v.ufl_index_dimensions
    total_shape = shape + index_shape
    return (shape, index_shape, total_shape)


def total_shape(v):
    """Compute the total shape of an expr."""
    sh, ish, tsh = compute_all_shapes(v)
    return tsh
