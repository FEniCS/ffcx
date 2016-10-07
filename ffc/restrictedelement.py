# -*- coding: utf-8 -*-
# Copyright (C) 2010 Anders Logg
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
# Modified by Marie Rognes, 2010.

import numpy
from ufl.utils.sorting import sorted_by_key
from ffc.log import error


class RestrictedElement:

    "Create a restriction of a given FIAT element."

    def __init__(self, element, indices, restriction_domain):
        if len(indices) == 0:
            error("No point in creating empty RestrictedElement.")

        self._element = element
        self._indices = indices
        self._entity_dofs = _extract_entity_dofs(element, indices)
        self._restriction_domain = restriction_domain

    def space_dimension(self):
        return len(self._indices)

    def value_shape(self):
        return self._element.value_shape()

    def degree(self):
        return self._element.degree()

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
        for (dtuple, values) in sorted_by_key(result):
            extracted[dtuple] = numpy.array([values[i] for i in self._indices])
        return extracted

    # Used in evaluate_basis:
    def get_coeffs(self):
        coefficients = self._element.get_coeffs()
        return numpy.array([coefficients[i] for i in self._indices])

    def dmats(self):
        return self._element.dmats()

    def get_num_members(self, arg):
        return self._element.get_num_members(arg)

    def restriction_domain(self):
        return self._restriction_domain


def _extract_entity_dofs(element, indices):
    # FIXME: Readability counts
    entity_dofs = element.entity_dofs()
    dofs = {}
    for (dim, entities) in sorted_by_key(entity_dofs):
        dofs[dim] = {}
        for (entity, all_dofs) in sorted_by_key(entities):
            dofs[dim][entity] = []
            for index in all_dofs:
                if index in indices:
                    # print "index = ", index
                    i = indices.index(index)
                    dofs[dim][entity] += [i]
    return dofs
