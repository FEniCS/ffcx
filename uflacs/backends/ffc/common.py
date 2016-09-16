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

"""FFC specific utilities."""

from six.moves import xrange as range


def ufc_restriction_postfix(restriction):
    # TODO: Get restriction postfix from somewhere central
    if restriction == "+":
        res = "_0"
    elif restriction == "-":
        res = "_1"
    else:
        res = ""
    return res


def ufc_restriction_offset(restriction, length):
    # TODO: Get restriction postfix from somewhere central
    if restriction == "-":
        return length
    else:
        return 0


# FIXME: Do something like this for shared symbol naming?
class FFCBackendSymbols(object):
    def __init__(self, language, coefficient_numbering):
        self.L = language
        self.S = self.L.Symbol
        self.coefficient_numbering = coefficient_numbering

        # Used for padding variable names based on restriction
        self.restriction_postfix = { r: ufc_restriction_postfix(r)
                                     for r in ("+", "-", None) }

    # FIXME: Used in access: weights, points, ia, A, w, x, J

    def entity(self, entitytype, restriction):
        "Entity index."
        return self.S(format_entity_name(entitytype, restriction))

    def x(self, quadloop):
        "Physical coordinates."
        return self.S("x" + str(quadloop))

    def xi(self, quadloop):
        "Reference cell coordinates."
        return self.S("xi" + str(quadloop))

    def quadrature_loop_index(self):
        "Reusing a single index name for all quadrature loops, assumed not to be nested."
        # If we want to use num_points-specific names for any symbols, this need num_points as well (or some other scope id).
        #return self.S("iq%d" % (num_points,))
        return self.S("iq")

    def coefficient_dof_sum_index(self):
        "Reusing a single index name for all coefficient dof*basis sums, assumed to always be the innermost loop."
        return self.S("ic")

    def coefficient_value_access(self, coefficient):  # TODO: Currently not used?
        c = self.coefficient_numbering[coefficient] # coefficient.count()
        # If we want to use num_points-specific names for any symbols, this need num_points as well (or some other scope id).
        #return self.S("w%d_%d" % (c, num_points))
        return self.S("w%d" % c)

    def coefficient_dof_access(self, coefficient, dof_number):
        # TODO: Add domain_number = self.ir["domain_numbering"][coefficient.ufl_domain().domain_key()]
        # TODO: Flatten dofs array and use CRSArray lookup table?
        # TODO: Apply integral specific renumbering?
        c = self.coefficient_numbering[coefficient] # coefficient.count()
        w = self.S("w")
        return w[c, dof_number]

    def domain_dof_access(self, dof, component, gdim, num_scalar_dofs, restriction, interleaved_components):
        # TODO: Add domain number as argument here, and {domain_offset} to array indexing:
        # domain_offset = self.ir["domain_offsets"][domain_number]
        vc = self.S("coordinate_dofs" + ufc_restriction_postfix(restriction))
        if interleaved_components:
            return vc[gdim*dof + component]
        else:
            return vc[num_scalar_dofs*component + dof]

    def domain_dofs_access(self, gdim, num_scalar_dofs, restriction, interleaved_components):
        # TODO: Add domain number as argument here, and {domain_offset} to array indexing:
        # domain_offset = self.ir["domain_offsets"][domain_number]
        return [self.domain_dof_access(dof, component, gdim, num_scalar_dofs, restriction, interleaved_components)
                for component in range(gdim)
                for dof in range(num_scalar_dofs)]


# TODO: This is not used much anymore, integrate in backend class, and use L.Symbol
class Names:

    def __init__(self):
        # Topology argument names
        self.vertex = "vertex"
        self.facet = "facet"

        # Geometry names
        self.coordinate_dofs = "coordinate_dofs"
        self.xi = "xi"
        self.x = "x"
        self.J = "J"
        self.K = "K"
        self.detJ = "detJ"
        self.det = "det"

        # Quadrature rule
        self.points = "points"
        self.weights = "weights"

        # Quadrature temps
        self.qw = "qw"
        self.D = "D"

        # (Base)name for intermediate registers
        self.s = "s"

        # Element tensor
        self.A = "A"

        # Coefficient dofs array
        self.w = "w"

        # Basenames for function components
        self.wbase = "w"
        self.vbase = "v"
        self.dwbase = "dw"
        self.dvbase = "dv"

        # Loop indices
        self.iq = "iq"   # Quadrature loop
        self.ic = "ic"   # Coefficient accumulation loop
        self.ia = "ia"   # Argument dof loop
        self.ild = "ild"  # Local derivative accumulation loop

        # Used for padding variable names based on restriction
        self.restriction_postfix = { r: ufc_restriction_postfix(r)
                                     for r in ("+", "-", None) }

names = Names()


def format_entity_name(entitytype, r):
    if entitytype == "cell":
        entity = "0"  # None # TODO: Keep 3D tables and use entity 0 for cells or make tables 2D and use None?
    elif entitytype == "facet":
        entity = names.facet + names.restriction_postfix[r]
    elif entitytype == "vertex":
        entity = names.vertex
    return entity


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
    return names.restriction_postfix[mt.restriction].replace("_", "_r")


def format_mt_name(basename, mt):
    access = "{basename}{avg}{res}{der}{comp}".format(basename=basename,
                                                      avg=format_mt_avg(mt),
                                                      res=format_mt_res(mt),
                                                      der=format_mt_der(mt),
                                                      comp=format_mt_comp(mt))
    return access
