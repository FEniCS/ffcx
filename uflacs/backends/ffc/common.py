# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""FFC/UFC specific symbol naming."""


from ffc.log import error


physical_quadrature_integral_types = ("custom", "cutcell", "interface", "overlap")


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
    for entity_dim in range(tdim+1):
        # n = dofs per cell entity
        if entity_dim == 0:
            n = 1
        elif entity_dim == 1:
            n = degree - 1
        elif entity_dim == 2:
            n = (degree - 2)*(degree - 1) // 2
        elif entity_dim == 3:
            n = (degree - 3)*(degree - 2)*(degree - 1) // 6
        else:
            error("Entity dimension out of range")
        # Accumulate
        num_entities = num_cell_entities[cellname][entity_dim]
        d += num_entities * n
    return d


# TODO: Get restriction postfix from somewhere central
def ufc_restriction_postfix(restriction):
    if restriction == "+":
        res = "_0"
    elif restriction == "-":
        res = "_1"
    else:
        res = ""
    return res


# TODO: Get restriction postfix from somewhere central
def ufc_restriction_offset(restriction, length):
    if restriction == "-":
        return length
    else:
        return 0


def format_mt_der(mt):
    # Expecting only local derivatives here
    assert not mt.global_derivatives
    # Add derivatives to name
    if mt.local_derivatives:
        der = "_d{0}".format(''.join(map(str, mt.local_derivatives)))
    else:
        der = ""
    return der


def format_mt_comp(mt):
    # Add flattened component to name (TODO: this should be the local component?)
    if mt.component:
        comp = "_c{0}".format(mt.flat_component)
    else:
        comp = ""
    return comp


def format_mt_avg(mt):
    # Add averaged state to name
    if mt.averaged:
        avg = "_a{0}".format(mt.averaged)
    else:
        avg = ""
    return avg


def format_mt_res(mt):
    return ufc_restriction_postfix(mt.restriction).replace("_", "_r")


def format_mt_name(basename, mt):
    access = "{basename}{avg}{res}{der}{comp}".format(basename=basename,
                                                      avg=format_mt_avg(mt),
                                                      res=format_mt_res(mt),
                                                      der=format_mt_der(mt),
                                                      comp=format_mt_comp(mt))
    return access


class FFCBackendSymbols(object):
    """FFC specific symbol definitions. Provides non-ufl symbols."""
    def __init__(self, language, coefficient_numbering):
        self.L = language
        self.S = self.L.Symbol
        self.coefficient_numbering = coefficient_numbering

        # Used for padding variable names based on restriction
        self.restriction_postfix = { r: ufc_restriction_postfix(r)
                                     for r in ("+", "-", None) }

    def element_tensor(self):
        "Symbol for the element tensor itself."
        return self.S("A")

    def entity(self, entitytype, restriction):
        "Entity index for lookup in element tables."
        if entitytype == "cell":
            # Always 0 for cells (even with restriction)
            return self.L.LiteralInt(0)
        elif entitytype == "facet":
            return self.S("facet" + ufc_restriction_postfix(restriction))
        elif entitytype == "vertex":
            return self.S("vertex")
        else:
            error("Unknown entitytype {}".format(entitytype))

    def cell_orientation_argument(self, restriction):
        "Cell orientation argument in ufc. Not same as cell orientation in generated code."
        return self.S("cell_orientation" + ufc_restriction_postfix(restriction))

    def cell_orientation_internal(self, restriction):
        "Internal value for cell orientation in generated code."
        return self.S("co" + ufc_restriction_postfix(restriction))

    def argument_loop_index(self, iarg):
        "Loop index for argument #iarg."
        return self.S("ia%d" % (iarg,))

    def quadrature_loop_index(self, num_points):
        """Reusing a single index name for all quadrature loops,
        assumed not to be nested."""
        if num_points == 1:
            return self.L.LiteralInt(0)
        elif num_points is None:
            return self.S("iq")
        else:
            return self.S("iq" + str(num_points))

    def coefficient_dof_sum_index(self):
        """Reusing a single index name for all coefficient dof*basis sums,
        assumed to always be the innermost loop."""
        return self.S("ic")

    def weights_array(self, num_points):
        return self.S("weights%d" % (num_points,))

    def points_array(self, num_points):
        # Note: Points array refers to points on the integration cell
        return self.S("points%d" % (num_points,))

    def physical_quadrature_points_array(self):
        return self.S("quadrature_points")

    def x_component(self, mt):
        "Physical coordinate component."
        return self.S(format_mt_name("x", mt))

    def J_component(self, mt):
        "Jacobian component."
        return self.S(format_mt_name("J", mt))

    def domain_dof_access(self, dof, component, gdim, num_scalar_dofs,
                          restriction, interleaved_components):
        # TODO: Add domain number?
        vc = self.S("coordinate_dofs" + ufc_restriction_postfix(restriction))
        if interleaved_components:
            return vc[gdim*dof + component]
        else:
            return vc[num_scalar_dofs*component + dof]

    def domain_dofs_access(self, gdim, num_scalar_dofs, restriction,
                           interleaved_components):
        # TODO: Add domain number?
        return [self.domain_dof_access(dof, component, gdim, num_scalar_dofs,
                                       restriction, interleaved_components)
                for component in range(gdim)
                for dof in range(num_scalar_dofs)]

    def coefficient_dof_access(self, coefficient, dof_number):
        # TODO: Add domain number?
        c = self.coefficient_numbering[coefficient]
        w = self.S("w")
        return w[c, dof_number]

    def coefficient_value(self, mt, num_points):
        "Symbol for variable holding value or derivative component of coefficient."
        c = self.coefficient_numbering[mt.terminal]
        return self.S(format_mt_name("w%d" % (c,), mt))
        # TODO: Should we include num_points here? Not sure if there is a need.
        #return self.S(format_mt_name("w%d_%d" % (c, num_points), mt))
