__author__ = "Anders Logg (logg@simula.no) and friends"
__copyright__ = "Copyright (C) 2005-2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006-2009
# Modified by Marie E. Rognes (meg@simula.no) 2007--2010
# Modified by Kristian B. Oelgaard 2009

# Last changed: 2010-01-14

# Python modules
import numpy

# UFL modules
from ufl import MixedElement as UFLMixedElement

# FFC modules
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

    def tabulate(self, order, points):
        """Tabulate values on mixed element by appropriately reordering
        the tabulated values for the sub elements."""

        # Special case: only one element
        if len(self._elements) == 1:
            return elements[0].tabulate(order, points)

        # Iterate over sub elements and build mixed table from element tables
        mixed_table = []
        offset = 0
        for element in self._elements:

            # Tabulate element
            element_table = element.tabulate(order, points)

            # Iterate over the components corresponding to the current element
            value_rank = len(element.value_shape())
            if value_rank == 0:
                component_table = _compute_component_table(element_table, offset, self.space_dimension())
                mixed_table.append(component_table)
            else:
                for i in range(element.value_dimension(0)):
                    component_table = _compute_component_table(element_table[i], offset, self.space_dimension())
                    mixed_table.append(component_table)

            # Add to offset, the number of the first basis function for the current element
            offset += element.space_dimension()

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

def _compute_component_table(table, offset, space_dimension):
    "Compute subtable for given component."

    component_table = {}

    # Iterate over derivative tuples
    for dtuple in table:

        # Get subtable for derivative tuple
        element_subtable = table[dtuple]

        # Initialize mixed subtable with zeros
        num_points = numpy.shape(element_subtable)[1]
        mixed_subtable = numpy.zeros((space_dimension, num_points), dtype = numpy.float)

        # Iterate over element basis functions and fill in non-zero values
        for i in range(len(element_subtable)):
            mixed_subtable[offset + i] = element_subtable[i]

        # Add to dictionary
        component_table[dtuple] = mixed_subtable

    return component_table
