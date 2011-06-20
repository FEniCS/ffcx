# Copyright (C) 2005-2010 Anders Logg
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
# Modified by Garth N. Wells, 2006-2009
# Modified by Marie E. Rognes, 2007-2010
# Modified by Kristian B. Oelgaard, 2010
#
# Last changed: 2010-01-30

# Python modules
import numpy

# FFC modules
from ffc.log import error

class MixedElement:
    "Create a FFC mixed element from a list of FFC/FIAT elements."

    def __init__(self, elements):
        self._elements = elements
        self._entity_dofs = _combine_entity_dofs(self._elements)

    def elements(self):
        return self._elements

    def space_dimension(self):
        return sum(e.space_dimension() for e in self._elements)

    def value_shape(self):
        # FIXME: value_shape for tensor elements in mixed elements not
        # well-defined
        return (sum(sum(e.value_shape()) or 1 for e in self._elements),)

    def entity_dofs(self):
        return self._entity_dofs

    def mapping(self):
        return [m for e in self._elements for m in e.mapping()]

    def dual_basis(self):
        return [L for e in self._elements for L in e.dual_basis()]

    def num_components(self):
        return sum(_num_components(e) for e in self._elements)

    def tabulate(self, order, points):
        """
        Tabulate values on mixed element by appropriately reordering
        the tabulated values for the nested elements.

        The table of values is organized as follows:

          D^a v_i[j](x_k) = table[a][i][j][k]

        where a is a multi-index (tuple) representing the derivative.
        For example, a = (1, 1) represents d/dx d/dy.
        """

        # Special case: only one element
        # NOTE: KBO: Why do we need this special case? (FFC Bug #798578)
        #       When calling tabulate() on a MixedElement one should be able to
        #       rely on getting data back which is ordered like a mixed element
        #       irrespective of the number of elements?
#        if len(self._elements) == 1:
#            return self._elements[0].tabulate(order, points)

        # Zeros for insertion into mixed table
        table_shape = (self.space_dimension(), self.num_components(), len(points))

        # Iterate over elements and fill in non-zero values
        irange = (0, 0)
        crange = (0, 0)
        mixed_table = {}
        for element in self._elements:

            # Tabulate element
            table = element.tabulate(order, points)

            # Compute range for insertion into mixed table
            irange = (irange[1], irange[1] + element.space_dimension())
            crange = (crange[1], crange[1] + _num_components(element))

            # Insert table into mixed table
            for dtuple in table.keys():

                # Insert zeros if necessary (should only happen first time)
                if not dtuple in mixed_table:
                    # NOTE: It is super important to create a new numpy.zeros
                    # instance to avoid manipulating a numpy reference in case
                    # it is created outside the loop.
                    mixed_table[dtuple] = numpy.zeros(table_shape)

                # Insert non-zero values
                if (crange[1] - crange[0]) > 1:
                    mixed_table[dtuple][irange[0]:irange[1], crange[0]:crange[1]] = table[dtuple]
                else:
                    mixed_table[dtuple][irange[0]:irange[1], crange[0]] = table[dtuple]

        return mixed_table

#--- Utility functions ---

def _combine_entity_dofs(elements):
    """
    Combine the entity_dofs from a list of elements into a combined
    entity_dof containing the information for all the elements.
    """

    # Return {} if no elements
    if not elements:
        return {}

    # Initialize entity_dofs dictionary
    entity_dofs = dict((key, {}) for key in elements[0].entity_dofs())
    for dim in elements[0].entity_dofs():
        for entity in elements[0].entity_dofs()[dim]:
            entity_dofs[dim][entity] = []

    offset = 0

    # Insert dofs from each element into the mixed entity_dof.
    for e in elements:
        dofs = e.entity_dofs()
        for dim in dofs:
            for entity in dofs[dim]:

                # Must take offset into account
                shifted_dofs = [v + offset for v in dofs[dim][entity]]

                # Insert dofs from this element into the entity_dofs
                entity_dofs[dim][entity] += shifted_dofs

        # Adjust offset
        offset += e.space_dimension()
    return entity_dofs

def _num_components(element):
    "Compute number of components for element."
    num_components = 0
    value_shape = element.value_shape()
    if len(value_shape) == 0:
        num_components += 1
    elif len(value_shape) == 1:
        num_components += value_shape[0]
    else:
        error("Cannot handle tensor-valued elements.")
    return num_components
