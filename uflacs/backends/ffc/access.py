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

from ufl.permutation import build_component_numbering
from ufl.corealg.multifunction import MultiFunction
from ufl.checks import is_cellwise_constant
from ffc.log import error
from ffc.log import ffc_assert

from uflacs.backends.ffc.common import FFCBackendSymbols

# FIXME: Move these to FFCBackendSymbols
from uflacs.backends.ffc.common import (names,
                                        format_entity_name,
                                        format_mt_name)


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

        # Need this for custom integrals
        #classname = make_classname(prefix, "finite_element", ir["element_numbers"][ufl_element])

        coefficient_numbering = self.ir["uflacs"]["coefficient_numbering"]
        self.symbols = FFCBackendSymbols(self.language, coefficient_numbering)

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

    def physical_points_array_name(self):
        return names.points

    def quadrature_loop_index(self):
        return self.symbols.quadrature_loop_index()

    def argument_loop_index(self, iarg):
        L = self.language
        return L.Symbol("{name}{num}".format(name=names.ia, num=iarg))

    def element_tensor_name(self):
        return names.A

    def element_tensor_entry(self, indices, shape):
        L = self.language
        flat_index = L.flattened_indices(indices, shape)
        A = L.Symbol(self.element_tensor_name())
        return A[flat_index]


    # === Rules for all modified terminal types ===

    def expr(self, e, mt, tabledata, num_points):
        error("Missing handler for type {0}.".format(e._ufl_class_.__name__))


    # === Rules for literal constants ===

    def zero(self, e, mt, tabledata, num_points):
        # We shouldn't have derivatives of constants left at this point
        assert not (mt.global_derivatives or mt.local_derivatives)
        # NB! UFL doesn't retain float/int type information for zeros...
        L = self.language
        return L.LiteralFloat(0.0)

    def int_value(self, e, mt, tabledata, num_points):
        # We shouldn't have derivatives of constants left at this point
        assert not (mt.global_derivatives or mt.local_derivatives)
        L = self.language
        return L.LiteralInt(int(e))

    def float_value(self, e, mt, tabledata, num_points):
        # We shouldn't have derivatives of constants left at this point
        assert not (mt.global_derivatives or mt.local_derivatives)
        L = self.language
        return L.LiteralFloat(float(e))

    def argument(self, e, mt, tabledata, num_points):
        L = self.language
        # Expecting only local derivatives and values here
        assert not mt.global_derivatives
        # assert mt.global_component is None

        # No need to store basis function value in its own variable, just get table value directly
        uname, begin, end = tabledata
        uname = L.Symbol(uname)

        entity = format_entity_name(self.ir["entitytype"], mt.restriction)
        entity = L.Symbol(entity)

        iq = self.quadrature_loop_index()
        idof = self.argument_loop_index(mt.terminal.number())

        return uname[entity][iq][idof - begin]

    def coefficient(self, e, mt, tabledata, num_points):
        t = mt.terminal
        if is_cellwise_constant(t):
            access = self._constant_coefficient(e, mt, tabledata)
        else:
            access = self._varying_coefficient(e, mt, tabledata)
        return access

    def _constant_coefficient(self, e, mt, tabledata):
        # Map component to flat index
        vi2si, si2vi = build_component_numbering(mt.terminal.ufl_shape,
                                                 mt.terminal.ufl_element().symmetry())
        ffc_assert(mt.flat_component == vi2si[mt.component], "Incompatible component flattening!")
        num_flat_components = len(si2vi)

        # Offset index if on second cell in interior facet integral
        offset = ufc_restriction_offset(mt.restriction, num_flat_components)
        if offset:
            idof = mt.flat_component + offset

        # Return direct reference to dof array
        return self.symbols.coefficient_dof_access(mt.terminal, idof)

    def _varying_coefficient(self, e, mt, tabledata):
        # Format base coefficient (derivative) name
        L = self.language
        coefficient_numbering = self.ir["uflacs"]["coefficient_numbering"]
        c = coefficient_numbering[mt.terminal] # mt.terminal.count()
        basename = "{name}{count}".format(name=names.w, count=c)
        return L.Symbol(format_mt_name(basename, mt))

    def quadrature_weight(self, e, mt, tabledata, num_points):
        L = self.language
        weight = self.weights_array_name(num_points)
        weight = L.Symbol(weight)
        iq = self.quadrature_loop_index()
        return weight[iq]

    def spatial_coordinate(self, e, mt, tabledata, num_points):
        L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of SpatialCoordinates.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of SpatialCoordinates.")
        #ffc_assert(not mt.restriction, "Not expecting restriction of SpatialCoordinates.")
        ffc_assert(not mt.averaged, "Not expecting average of SpatialCoordinates.")

        if self.physical_coordinates_known:
            # In a context where the physical coordinates are available in existing variables.
            x = self.physical_points_array_name()
            x = L.Symbol(x)
            iq = self.quadrature_loop_index()
            gdim, = mt.terminal.ufl_shape
            return x[iq * gdim + mt.flat_component]
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
            X = self.points_array_name(num_points)
            X = L.Symbol(X)
            iq = self.quadrature_loop_index()
            tdim, = mt.terminal.ufl_shape
            return X[iq * tdim + mt.flat_component]

    def jacobian(self, e, mt, tabledata, num_points):
        L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of Jacobian.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of Jacobian.")
        ffc_assert(not mt.averaged, "Not expecting average of Jacobian.")

        return L.Symbol(format_mt_name(names.J, mt))

    def reference_cell_volume(self, e, mt, tabledata, access):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            return L.Symbol("{0}_reference_cell_volume".format(cellname))
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def reference_facet_volume(self, e, mt, tabledata, access):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            return L.Symbol("{0}_reference_facet_volume".format(cellname))
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def reference_normal(self, e, mt, tabledata, access):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol("{0}_reference_facet_normals".format(cellname))
            facet = L.Symbol(format_entity_name("facet", mt.restriction))
            return table[facet][mt.component[0]]
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def cell_facet_jacobian(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol("{0}_reference_facet_jacobian".format(cellname))
            facet = L.Symbol(format_entity_name("facet", mt.restriction))
            return table[facet][mt.component[0]][mt.component[1]]
        elif cellname == "interval":
            error("The reference facet jacobian doesn't make sense for interval cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def cell_edge_vectors(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol("{0}_reference_edge_vectors".format(cellname))
            return table[mt.component[0]][mt.component[1]]
        elif cellname == "interval":
            error("The reference cell edge vectors doesn't make sense for interval cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def facet_edge_vectors(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("tetrahedron", "hexahedron"):
            table = L.Symbol("{0}_reference_edge_vectors".format(cellname))
            facet = L.Symbol(format_entity_name("facet", mt.restriction))
            return table[facet][mt.component[0]][mt.component[1]]
        elif cellname in ("interval", "triangle", "quadrilateral"):
            error("The reference cell facet edge vectors doesn't make sense for interval or triangle cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))

    def cell_orientation(self, e, mt, tabledata, num_points):
        L = self.language
        # Error if not in manifold case:
        gdim = mt.terminal.ufl_domain().geometric_dimension()
        tdim = mt.terminal.ufl_domain().topological_dimension()
        assert gdim > tdim
        return L.Symbol("co")

    def facet_orientation(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname not in ("interval", "triangle", "tetrahedron"):
            error("Unhandled cell types {0}.".format(cellname))

        table = L.Symbol("{0}_facet_orientations".format(cellname))
        facet = L.Symbol(format_entity_name("facet", mt.restriction))
        return table[facet]

    def _expect_symbolic_lowering(self, e, mt, tabledata, num_points):
        error("Expecting {0} to be replaced in symbolic preprocessing.".format(type(e)))
    facet_normal = _expect_symbolic_lowering
    cell_normal = _expect_symbolic_lowering
    jacobian_inverse = _expect_symbolic_lowering
    jacobian_determinant = _expect_symbolic_lowering
    facet_jacobian = _expect_symbolic_lowering
    facet_jacobian_inverse = _expect_symbolic_lowering
    facet_jacobian_determinant = _expect_symbolic_lowering
