# -*- coding: utf-8 -*-
# Copyright (C) 2011-2016 Martin Sandve Alnæs
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

"""FFC/UFC specific variable access."""

from ufl.corealg.multifunction import MultiFunction
from ufl.permutation import build_component_numbering
from ufl.measure import custom_integral_types

from ffc.log import error, warning, debug

from ffc.uflacs.backends.ffc.symbols import FFCBackendSymbols


class FFCBackendAccess(MultiFunction):
    """FFC specific cpp formatter class."""

    def __init__(self, ir, language, symbols, parameters):
        MultiFunction.__init__(self)

        # Store ir and parameters
        self.entitytype = ir["entitytype"]
        self.integral_type = ir["integral_type"]
        self.language = language
        self.symbols = symbols
        self.parameters = parameters


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


    #def quadrature_weight(self, e, mt, tabledata, num_points):
    #    "Quadrature weights are precomputed and need no code."
    #    return []


    def coefficient(self, e, mt, tabledata, num_points):
        ttype = tabledata.ttype

        assert ttype != "zeros"

        begin, end = tabledata.dofrange

        if ttype == "ones" and (end - begin) == 1:
            # f = 1.0 * f_{begin}, just return direct reference to dof array at dof begin
            # (if mt is restricted, begin contains cell offset)
            idof = begin
            return self.symbols.coefficient_dof_access(mt.terminal, idof)
        elif ttype == "quadrature":
            # Dofmap should be contiguous in this case
            assert len(tabledata.dofmap) == end - begin
            # f(x_q) = sum_i f_i * delta_iq = f_q, just return direct
            # reference to dof array at quadrature point index + begin
            iq = self.symbols.quadrature_loop_index()
            idof = begin + iq
            return self.symbols.coefficient_dof_access(mt.terminal, idof)
        else:
            # Return symbol, see definitions for computation 
            return self.symbols.coefficient_value(mt)  #, num_points)


    def spatial_coordinate(self, e, mt, tabledata, num_points):
        #L = self.language
        if mt.global_derivatives:
            error("Not expecting global derivatives of SpatialCoordinate.")
        if mt.averaged:
            error("Not expecting average of SpatialCoordinates.")

        if self.integral_type in custom_integral_types:
            if mt.local_derivatives:
                error("FIXME: Jacobian in custom integrals is not implemented.")

            # Access predefined quadrature points table
            x = self.symbols.custom_points_table()
            iq = self.symbols.quadrature_loop_index()
            gdim, = mt.terminal.ufl_shape
            if gdim == 1:
                index = iq
            else:
                index = iq * gdim + mt.flat_component
            return x[index]
        else:
            # Physical coordinates are computed by code generated in definitions
            return self.symbols.x_component(mt)


    def cell_coordinate(self, e, mt, tabledata, num_points):
        #L = self.language
        if mt.global_derivatives:
            error("Not expecting derivatives of CellCoordinate.")
        if mt.local_derivatives:
            error("Not expecting derivatives of CellCoordinate.")
        if mt.averaged:
            error("Not expecting average of CellCoordinate.")

        if self.integral_type == "cell" and not mt.restriction:
            # Access predefined quadrature points table
            X = self.symbols.points_table(num_points)
            tdim, = mt.terminal.ufl_shape
            iq = self.symbols.quadrature_loop_index()
            if num_points == 1:
                index = mt.flat_component
            elif tdim == 1:
                index = iq
            else:
                index = iq * tdim + mt.flat_component
            return X[index]
        else:
            # X should be computed from x or Xf symbolically instead of getting here
            error("Expecting reference cell coordinate to be symbolically rewritten.")


    def facet_coordinate(self, e, mt, tabledata, num_points):
        L = self.language
        if mt.global_derivatives:
            error("Not expecting derivatives of FacetCoordinate.")
        if mt.local_derivatives:
            error("Not expecting derivatives of FacetCoordinate.")
        if mt.averaged:
            error("Not expecting average of FacetCoordinate.")
        if mt.restriction:
            error("Not expecting restriction of FacetCoordinate.")

        if self.integral_type in ("interior_facet", "exterior_facet"):
            tdim, = mt.terminal.ufl_shape
            if tdim == 0:
                error("Vertices have no facet coordinates.")
            elif tdim == 1:
                # 0D vertex coordinate
                warning("Vertex coordinate is always 0, should get rid of this in ufl geometry lowering.")
                return L.LiteralFloat(0.0)
            Xf = self.points_table(num_points)
            iq = self.symbols.quadrature_loop_index()
            assert 0 <= mt.flat_component < (tdim-1)
            if num_points == 1:
                index = mt.flat_component
            elif tdim == 2:
                index = iq
            else:
                index = iq * (tdim - 1) + mt.flat_component
            return Xf[index]
        else:
            # Xf should be computed from X or x symbolically instead of getting here
            error("Expecting reference facet coordinate to be symbolically rewritten.")


    def jacobian(self, e, mt, tabledata, num_points):
        L = self.language
        if mt.global_derivatives:
            error("Not expecting global derivatives of Jacobian.")
        if mt.averaged:
            error("Not expecting average of Jacobian.")
        return self.symbols.J_component(mt)


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
            facet = self.symbols.entity("facet", mt.restriction)
            return table[facet][mt.component[0]]
        else:
            error("Unhandled cell types {0}.".format(cellname))


    def cell_facet_jacobian(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol("{0}_reference_facet_jacobian".format(cellname))
            facet = self.symbols.entity("facet", mt.restriction)
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
            facet = self.symbols.entity("facet", mt.restriction)
            return table[facet][mt.component[0]][mt.component[1]]
        elif cellname in ("interval", "triangle", "quadrilateral"):
            error("The reference cell facet edge vectors doesn't make sense for interval or triangle cell.")
        else:
            error("Unhandled cell types {0}.".format(cellname))


    def cell_orientation(self, e, mt, tabledata, num_points):
        # Error if not in manifold case:
        domain = mt.terminal.ufl_domain()
        assert domain.geometric_dimension() > domain.topological_dimension()
        return self.symbols.cell_orientation_internal(mt.restriction)


    def facet_orientation(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname not in ("interval", "triangle", "tetrahedron"):
            error("Unhandled cell types {0}.".format(cellname))

        table = L.Symbol("{0}_facet_orientations".format(cellname))
        facet = self.symbols.entity("facet", mt.restriction)
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
