# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFC/UFC specific symbol naming."""

import logging

logger = logging.getLogger(__name__)


# TODO: Move somewhere else
def num_coordinate_component_dofs(coordinate_element):
    """Get the number of dofs for a coordinate component for this degree.

    This is a local hack that works for Lagrange 1-3, better
    would be to get this passed by ffc from fiat through the ir.
    The table data is to messy to figure out a clean design for that quickly.
    """
    from ufl.cell import num_cell_entities
    degree = coordinate_element.degree()
    cell = coordinate_element.cell()
    tdim = cell.topological_dimension()
    cellname = cell.cellname()
    d = 0
    for entity_dim in range(tdim + 1):
        # n = dofs per cell entity
        if entity_dim == 0:
            n = 1
        elif entity_dim == 1:
            n = degree - 1
        elif entity_dim == 2:
            n = (degree - 2) * (degree - 1) // 2
        elif entity_dim == 3:
            n = (degree - 3) * (degree - 2) * (degree - 1) // 6
        else:
            raise RuntimeError("Entity dimension out of range")
        # Accumulate
        num_entities = num_cell_entities[cellname][entity_dim]
        d += num_entities * n
    return d


# TODO: Get restriction postfix from somewhere central
def ufc_restriction_offset(restriction, length):
    if restriction == "-":
        return length
    else:
        return 0
