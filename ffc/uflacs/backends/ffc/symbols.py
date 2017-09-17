# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
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


# TODO: Get restriction postfix from somewhere central
def ufc_restriction_postfix(restriction):
    if restriction == "+":
        res = "_0"
    elif restriction == "-":
        res = "_1"
    else:
        res = ""
    return res


def format_mt_name(basename, mt):
    "Format variable name for modified terminal."
    access = str(basename)

    # Add averaged state to name
    if mt.averaged:
        avg = "_a{0}".format(mt.averaged)
        access += avg

    # Format restriction
    res = ufc_restriction_postfix(mt.restriction).replace("_", "_r")
    access += res

    # Format local derivatives
    assert not mt.global_derivatives
    if mt.local_derivatives:
        der = "_d{0}".format(''.join(map(str, mt.local_derivatives)))
        access += der

    # Add flattened component to name
    if mt.component:
        comp = "_c{0}".format(mt.flat_component)
        access += comp

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

        # TODO: Make this configurable for easy experimentation with dolfin!
        # Coordinate dofs for each component are interleaved? Must match dolfin.
        # True = XYZXYZXYZXYZ, False = XXXXYYYYZZZZ
        self.interleaved_components = True


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
        indices = ["i", "j", "k", "l"]
        return self.S(indices[iarg])


    def coefficient_dof_sum_index(self):
        """Index for loops over coefficient dofs, assumed to never be used in two nested loops."""
        return self.S("ic")


    def quadrature_loop_index(self):
        """Reusing a single index name for all quadrature loops,
        assumed not to be nested."""
        return self.S("iq")


    def num_custom_quadrature_points(self):
        "Number of quadrature points, argument to custom integrals."
        return self.S("num_quadrature_points")


    def custom_quadrature_weights(self):
        "Quadrature weights including cell measure scaling, argument to custom integrals."
        return self.S("quadrature_weights")


    def custom_quadrature_points(self):
        "Physical quadrature points, argument to custom integrals."
        return self.S("quadrature_points")


    def custom_weights_table(self):
        "Table for chunk of custom quadrature weights (including cell measure scaling)."
        return self.S("weights_chunk")


    def custom_points_table(self):
        "Table for chunk of custom quadrature points (physical coordinates)."
        return self.S("points_chunk")


    def weights_table(self, num_points):
        "Table of quadrature weights."
        return self.S("weights%d" % (num_points,))


    def points_table(self, num_points):
        "Table of quadrature points (points on the reference integration entity)."
        return self.S("points%d" % (num_points,))


    def x_component(self, mt):
        "Physical coordinate component."
        return self.S(format_mt_name("x", mt))


    def X_component(self, mt):
        "Reference coordinate component."
        return self.S(format_mt_name("X", mt))


    def J_component(self, mt):
        "Jacobian component."
        # FIXME: Add domain number!
        return self.S(format_mt_name("J", mt))


    def domain_dof_access(self, dof, component, gdim, num_scalar_dofs, restriction):
        # FIXME: Add domain number or offset!
        vc = self.S("coordinate_dofs" + ufc_restriction_postfix(restriction))
        if self.interleaved_components:
            return vc[gdim*dof + component]
        else:
            return vc[num_scalar_dofs*component + dof]


    def domain_dofs_access(self, gdim, num_scalar_dofs, restriction):
        # FIXME: Add domain number or offset!
        return [self.domain_dof_access(dof, component, gdim, num_scalar_dofs, restriction)
                for component in range(gdim)
                for dof in range(num_scalar_dofs)]


    def coefficient_dof_access(self, coefficient, dof_number):
        # TODO: Add domain number?
        c = self.coefficient_numbering[coefficient]
        w = self.S("w")
        return w[c, dof_number]


    def coefficient_value(self, mt):
        "Symbol for variable holding value or derivative component of coefficient."
        c = self.coefficient_numbering[mt.terminal]
        return self.S(format_mt_name("w%d" % (c,), mt))


    def element_table(self, tabledata, entitytype, restriction):
        if tabledata.is_uniform:
            entity = 0
        else:
            entity = self.entity(entitytype, restriction)

        if tabledata.is_piecewise:
            iq = 0
        else:
            iq = self.quadrature_loop_index()

        # Return direct access to element table
        return self.S(tabledata.name)[entity][iq]
