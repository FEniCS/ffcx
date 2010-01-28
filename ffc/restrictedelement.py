__author__ = "Anders Logg (logg@simula.no) and friends"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-28

class RestrictedElement:
    "Create a restriction of a given element."

    def __init__(self, element):
        self._element = element

    def space_dimension(self):
        return self._element.space_dimension()

    def value_shape(self):
        return self._element.value_shape()

    def entity_dofs(self):
        return self._element.entity_dofs()

    def mapping(self):
        return self._element.mapping()

    def dual_basis(self):
        return self._element.dual_basis()

    def tabulate(self, order, points):
        return self._element.tabulate(order, points)
