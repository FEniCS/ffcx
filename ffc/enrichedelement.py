# Copyright (C) 2010 Marie E. Rognes
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
# First added:  2010-03-07
# Last changed: 2010-03-07

import numpy
from utils import pick_first
from mixedelement import _combine_entity_dofs, _num_components

class EnrichedElement:
    "Create the space spanned by a list of ffc elements."

    def __init__(self, elements):
        self._elements = elements
        self._entity_dofs = _combine_entity_dofs(elements)

    def elements(self):
        return self._elements

    def space_dimension(self):
        return sum(e.space_dimension() for e in self._elements)

    def value_shape(self):
        return pick_first([e.value_shape() for e in self._elements])

    def degree(self):
        return max(e.degree() for e in self._elements)

    def entity_dofs(self):
        return self._entity_dofs

    def mapping(self):
        return [m for e in self._elements for m in e.mapping()]

    def dual_basis(self):
        return [L for e in self._elements for L in e.dual_basis()]

    def tabulate(self, order, points):

        num_components = _num_components(self)
        table_shape = (self.space_dimension(), num_components, len(points))

        table = {}
        irange = (0, 0)
        for element in self._elements:

            etable = element.tabulate(order, points)
            irange = (irange[1], irange[1] + element.space_dimension())

            # Insert element table into table
            for dtuple in etable.keys():

                if not dtuple in table:
                    if num_components == 1:
                        table[dtuple] = numpy.zeros((self.space_dimension(), len(points)))
                    else:
                        table[dtuple] = numpy.zeros(table_shape)

                table[dtuple][irange[0]:irange[1]][:] = etable[dtuple]

        return table


class SpaceOfReals:

    def __init__(self, element):
        self._element = element
        self._entity_dofs = element.entity_dofs()

    def space_dimension(self):
        return 1

    def value_shape(self):
        return ()

    def degree(self):
        return 0

    def entity_dofs(self):
        return self._entity_dofs

    def mapping(self):
        return ["affine"]

    def dual_basis(self):
        return self._element.dual_basis()

    def tabulate(self, order, points):
        return self._element.tabulate(order, points)

    def get_coeffs(self):
        return self._element.get_coeffs()

    def dmats(self):
        return self._element.dmats()

    def get_num_members(self, arg):
        return self._element.get_num_members(arg)
