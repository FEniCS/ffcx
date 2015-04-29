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

"""FFC specific access formatting."""

import ufl
from ufl.permutation import build_component_numbering
from ufl.corealg.multifunction import MultiFunction

from ffc.log import error
from ffc.log import ffc_assert

# FIXME: These need language input:
from uflacs.backends.ffc.common import (names,
                                        format_entity_name,
                                        format_mt_name,
                                        generate_coefficient_dof_access)


class FFCAccessBackend(MultiFunction):
    """FFC specific cpp formatter class."""

    def __init__(self, ir, language, parameters):
        MultiFunction.__init__(self)

        # Store ir and parameters
        self.ir = ir
        self.language = language
        self.parameters = parameters

        # Configure definitions behaviour
        self.physical_coordinates_known = self.ir["integral_type"] == "quadrature"

    def get_includes(self):
        "Return include statements to insert at top of file."
        includes = []
        return includes


    # === Access to names of quantities not among the symbolic UFL types ===
    # FIXME: Move these out of the AccessBackend, maybe introduce a FFCBackendSymbols?
    #        A symbols class can contain generate*names from common.* as well.
    # FIXME: Use self.language.Symbol and/or self.language.ArrayAccess to wrap names.*:
    def weights_array_name(self, num_points):
        return "{0}{1}".format(names.weights, num_points)

    def points_array_name(self, num_points):
        return "{0}{1}".format(names.points, num_points)

    def quadrature_loop_index(self):  # (num_points):
        # If we want to use num_points-specific names for the quadrature loop index, definitions.py need num_points as well.
        # return "{0}{1}".format(names.iq, num_points)
        return names.iq

    def argument_loop_index(self, iarg):
        return "{name}{num}".format(name=names.ia, num=iarg)

    def element_tensor_name(self):
        return names.A

    def element_tensor_entry(self, indices, shape):
        L = self.language
        flat_index = L.flattened_indices(indices, shape)
        return L.ArrayAccess(names.A, flat_index)


    # === Rules for all modified terminal types ===

    def expr(self, e, mt, tabledata, num_points):
        error("Missing handler for type {0}.".format(e._ufl_class_.__name__))


    # === Rules for literal constants ===

    def zero(self, e, mt, tabledata, num_points):
        # We shouldn't have derivatives of constants left at this point
        assert not (mt.global_derivatives or mt.local_derivatives)
        # NB! UFL doesn't retain float/int type information for zeros...
        L = self.language
        return L.Literal(0.0)

    def int_value(self, e, mt, tabledata, num_points):
        # We shouldn't have derivatives of constants left at this point
        assert not (mt.global_derivatives or mt.local_derivatives)
        L = self.language
        return L.Literal(int(e))

    def float_value(self, e, mt, tabledata, num_points):
        # We shouldn't have derivatives of constants left at this point
        assert not (mt.global_derivatives or mt.local_derivatives)
        L = self.language
        return L.Literal(float(e))

    def argument(self, e, mt, tabledata, num_points):
        L = self.language
        # Expecting only local derivatives and values here
        assert not mt.global_derivatives
        # assert mt.global_component is None

        # No need to store basis function value in its own variable, just get table value directly
        uname, begin, end = tabledata
        entity = format_entity_name(self.ir["entitytype"], mt.restriction)

        iq = self.quadrature_loop_index()
        idof = self.argument_loop_index(mt.terminal.number())

        dof = L.Sub(idof, begin)
        return L.ArrayAccess(uname, (entity, iq, dof))

    def coefficient(self, e, mt, tabledata, num_points):
        t = mt.terminal
        if t.is_cellwise_constant():
            access = self._constant_coefficient(e, mt, tabledata)
        else:
            access = self._varying_coefficient(e, mt, tabledata)
        return access

    def _constant_coefficient(self, e, mt, tabledata):
        # Map component to flat index
        vi2si, si2vi = build_component_numbering(mt.terminal.ufl_shape,
                                                 mt.terminal.element().symmetry())
        num_flat_components = len(si2vi)
        ffc_assert(mt.flat_component == vi2si[mt.component], "Incompatible component flattening!")

        # Offset index if on second cell in interior facet integral
        # TODO: Get the notion that '-' is the second cell from a central definition?
        if mt.restriction == "-":
            idof = mt.flat_component + len(si2vi)
        else:
            idof = mt.flat_component

        # Return direct reference to dof array
        L = self.language
        return generate_coefficient_dof_access(L, mt.terminal, idof)

    def _varying_coefficient(self, e, mt, tabledata):
        # Format base coefficient (derivative) name
        basename = "{name}{count}".format(name=names.w, count=mt.terminal.count())
        L = self.language
        return L.Symbol(format_mt_name(basename, mt))

    def quadrature_weight(self, e, mt, tabledata, num_points):
        name = self.weights_array_name(num_points)
        ind = self.quadrature_loop_index()
        L = self.language
        return L.ArrayAccess(name, ind)

    def spatial_coordinate(self, e, mt, tabledata, num_points):
        L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of SpatialCoordinates.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of SpatialCoordinates.")
        # ffc_assert(not mt.restriction, "Not expecting restriction of SpatialCoordinates.")
        ffc_assert(not mt.averaged, "Not expecting average of SpatialCoordinates.")

        if self.physical_coordinates_known:
            # In a context where the physical coordinates are available in existing variables.
            name = self.points_array_name(num_points)
            ind = (self.quadrature_loop_index(), mt.flat_component)
            return L.ArrayAccess(name, ind)

        elif mt.terminal.domain().coordinates() is not None:
            # No special variable should exist in this case.
            error("Expecting spatial coordinate to be symbolically rewritten.")

        else:
            # In a context where physical coordinates are computed by code generated by us.
            return L.Symbol(format_mt_name(names.x, mt))

    def cell_coordinate(self, e, mt, tabledata, num_points):
        L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of CellCoordinates.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of CellCoordinates.")
        ffc_assert(not mt.averaged, "Not expecting average of CellCoordinates.")

        assert not mt.restriction  # FIXME: Not used!

        if self.physical_coordinates_known:
            # No special variable should exist in this case.
            error("Expecting reference coordinate to be symbolically rewritten.")
        else:
            name = self.points_array_name(num_points)
            ind = (self.quadrature_loop_index(), mt.flat_component)
            return L.ArrayAccess(name, ind)

    def jacobian(self, e, mt, tabledata, num_points):
        L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of Jacobian.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of Jacobian.")
        ffc_assert(not mt.averaged, "Not expecting average of Jacobian.")

        return L.Symbol(format_mt_name(names.J, mt))

    def cell_facet_jacobian(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.domain().cell().cellname()
        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            tablename = "{0}_reference_facet_jacobian".format(cellname)
            facet = format_entity_name("facet", mt.restriction)
            return L.ArrayAccess(tablename, (facet, mt.component[0], mt.component[1]))
        elif cellname == "interval":
            error("The reference facet jacobian doesn't make sense for interval cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def cell_edge_vectors(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.domain().cell().cellname()
        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            tablename = "{0}_reference_edge_vectors".format(cellname)
            return L.ArrayAccess(tablename, (mt.component[0], mt.component[1]))
        elif cellname == "interval":
            error("The reference cell edge vectors doesn't make sense for interval cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def facet_edge_vectors(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.domain().cell().cellname()
        if cellname in ("tetrahedron", "hexahedron"):
            tablename = "{0}_reference_edge_vectors".format(cellname)
            facet = format_entity_name("facet", mt.restriction)
            return L.ArrayAccess(tablename, (facet, mt.component[0], mt.component[1]))
        elif cellname in ("interval", "triangle", "quadrilateral"):
            error("The reference cell facet edge vectors doesn't make sense for interval or triangle cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def cell_orientation(self, e, mt, tabledata, num_points):
        L = self.language
        # Error if not in manifold case:
        gdim = mt.terminal.cell().geometric_dimension()
        tdim = mt.terminal.cell().topological_dimension()
        assert gdim > tdim
        return L.Symbol("cell_orientation")

    def facet_orientation(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.domain().cell().cellname()
        if cellname not in ("interval", "triangle", "tetrahedron"):
            error("Unhandled cell types {0}.".format(cellname))

        tablename = "{0}_facet_orientations".format(cellname)
        facet = format_entity_name("facet", mt.restriction)
        return L.ArrayAccess(tablename, (facet,))

    def _expect_symbolic_lowering(self, e, mt, tabledata, num_points):
        error("Expecting {0} to be replaced in symbolic preprocessing.".format(type(e)))
    facet_normal = _expect_symbolic_lowering
    cell_normal = _expect_symbolic_lowering
    jacobian_inverse = _expect_symbolic_lowering
    jacobian_determinant = _expect_symbolic_lowering
    facet_jacobian = _expect_symbolic_lowering
    facet_jacobian_inverse = _expect_symbolic_lowering
    facet_jacobian_determinant = _expect_symbolic_lowering
