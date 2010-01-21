__author__ = "Anders Logg (logg@simula.no) and friends"
__copyright__ = "Copyright (C) 2005-2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006-2009
# Modified by Marie E. Rognes (meg@simula.no) 2007--2010
# Modified by Kristian B. Oelgaard 2009

# Last changed: 2010-01-21

# Python modules
import numpy

# UFL modules
from ufl import MixedElement as UFLMixedElement

# FFC modules
from log import error
from fiatinterface import create_fiat_element

class MixedElement:
    """
    Create a FFC mixed element from a UFL mixed element.
    """

    def __init__(self, ufl_element):
        self._elements = _extract_elements(ufl_element)
        self._entity_dofs = _combine_entity_dofs(self._elements)
        self._value_shape = ufl_element.value_shape()

    def elements(self):
        return self._elements

    def space_dimension(self):
        return sum(e.space_dimension() for e in self._elements)

    def value_shape(self):
        return self._value_shape

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
        if len(self._elements) == 1:
            return elements[0].tabulate(order, points)

        # Zeros for insertion into mixed table
        table_shape = (self.space_dimension(), self.num_components(), len(points))
        zeros = numpy.zeros(table_shape)

        print "shape =", table_shape

        # Iterate over elements and fill in non-zero values
        irange = (0, 0)
        crange = (0, 0)
        mixed_table = {}
        for element in self._elements:

            print element
            print order

            # Tabulate element
            table = element.tabulate(order, points)

            # Compute range for insertion into mixed table
            irange = (irange[1], irange[1] + element.space_dimension())
            crange = (crange[1], crange[1] + _num_components(element))

            print "irange =", irange
            print "crange =", crange

            # Insert table into mixed table
            for dtuple in table.keys():

                print "shape  =", numpy.shape(table[dtuple])

                # Insert zeros if necessary (should only happen first time)
                if not dtuple in mixed_table:
                    mixed_table[dtuple] = zeros

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

def _extract_elements(ufl_element):
    "Recursively extract un-nested list of (component) elements."

    elements = []
    if isinstance(ufl_element, UFLMixedElement):
        for sub_element in ufl_element.sub_elements():
            elements += _extract_elements(sub_element)
        return elements

    elements += [create_fiat_element(ufl_element)]
    return elements

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
