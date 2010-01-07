__author__ = "Anders Logg (logg@simula.no) and friends"
__copyright__ = "Copyright (C) 2005-2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006-2009
# Modified by Marie E. Rognes (meg@simula.no) 2007--2010
# Modified by Kristian B. Oelgaard 2009

# Last changed: 2010-01-07

# FFC modules.
from fiatinterface import create_fiat_element

class MixedElement:
    """
    Create a FFC mixed element from a UFL mixed element
    """

    def __init__(self, ufl_element):
        self._elements = extract_elements(ufl_element)
        self._entity_dofs = combine_entity_dofs(self._elements)
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


# ---- Utility functions ----

def combine_entity_dofs(elements):
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

def extract_elements(ufl_element):
    # Q: Are UFL subelements nested?
    # Q: Vector elements pose some interesting issues...
    return [create_fiat_element(sub_element)
            for sub_element in ufl_element.sub_elements()]

