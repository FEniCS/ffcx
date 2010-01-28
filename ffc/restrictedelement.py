__author__ = "Anders Logg (logg@simula.no) and friends"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-29

import numpy

def _extract_entity_dofs(element, indices):

    # FIXME: Readability counts
    entity_dofs = element.entity_dofs()
    dofs = {}
    for (dim, entities) in entity_dofs.iteritems():
        dofs[dim] = {}
        for (entity, all_dofs) in entities.iteritems():
            dofs[dim][entity] = []
            for index in all_dofs:
                if index in indices:
                    i = indices.index(index)
                    dofs[dim][entity] += [i]
    return entity_dofs


class RestrictedElement:
    "Create a restriction of a given element."

    def __init__(self, element, indices):
        self._element = element
        self._indices = indices
        self._entity_dofs = _extract_entity_dofs(element, indices)

    def space_dimension(self):
        return len(self._indices)

    def value_shape(self):
        return self._element.value_shape()

    def entity_dofs(self):
        return self._entity_dofs

    def mapping(self):
        mappings = self._element.mapping()
        return [mappings[i] for i in self._indices]

    def dual_basis(self):
        dual = self._element.dual_basis()
        return [dual[i] for i in self._indices]

    def tabulate(self, order, points):
        result = self._element.tabulate(order, points)
        extracted = {}
        for (dtuple, values) in result.iteritems():
            extracted[dtuple] = numpy.array([values[i] for i in self._indices])
        return extracted
