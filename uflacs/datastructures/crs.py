# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""Compressed row storage 'matrix' (actually just a non-rectangular 2d array)."""

from six.moves import xrange as range
import numpy


class CRS(object):

    """A simple compressed row storage matrix.

    This CRS variant doesn't have a sparsity pattern,
    as each row is simply a dense vector.
    """

    def __init__(self, row_capacity, element_capacity, dtype):
        self.row_offsets = numpy.zeros(row_capacity + 1, dtype=int)
        self.data = numpy.zeros(element_capacity, dtype=dtype)
        self.num_rows = 0

    def push_row(self, elements):
        a = self.row_offsets[self.num_rows]
        b = a + len(elements)
        self.data[a:b] = elements
        self.num_rows += 1
        self.row_offsets[self.num_rows] = b

    @property
    def num_elements(self):
        return self.row_offsets[self.num_rows]

    def __getitem__(self, row):
        if row < 0 or row >= self.num_rows:
            raise IndexError("Row number out of range!")
        a = self.row_offsets[row]
        b = self.row_offsets[row + 1]
        return self.data[a:b]

    def __len__(self):
        return self.num_rows

    def __str__(self):
        return "[%s]" % (', '.join(str(row) for row in self),)


def list_to_crs(elements):
    "Construct a diagonal CRS matrix from a list of elements of the same type."
    n = len(elements)
    crs = CRS(n, n, type(elements[0]))
    for element in elements:
        crs.push_row((element,))
    return crs


def rows_dict_to_crs(rows, num_rows, num_elements, dtype):
    "Construct a CRS matrix from a dict mapping row index to row elements list."
    crs = CRS(num_rows, num_elements, dtype)
    for i in range(num_rows):
        crs.push_row(rows.get(i, ()))
    return crs


def rows_to_crs(rows, num_rows, num_elements, dtype):
    "Construct a CRS matrix from a list of row element lists."
    crs = CRS(num_rows, num_elements, dtype)
    for row in rows:
        crs.push_row(row)
    return crs
