__author__ = "Marie E. Rognes (meg@simula.no)"
__date__ = "2010-03-07"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from utils import pick_first
from mixedelement import _combine_entity_dofs

class ElementUnion:
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
