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

"""FFC/UFC specific variable access."""

from ufl.corealg.multifunction import MultiFunction
from ufl.permutation import build_component_numbering

from ffc.log import error, warning, ffc_assert

from uflacs.backends.ffc.common import FFCBackendSymbols
from uflacs.backends.ffc.common import physical_quadrature_integral_types


class FFCBackendAccess(MultiFunction):
    """FFC specific cpp formatter class."""

    def __init__(self, ir, language, symbols, parameters):
        MultiFunction.__init__(self)

        # Store ir and parameters
        self.ir = ir
        self.language = language
        self.symbols = symbols
        self.parameters = parameters

        # Need this for custom integrals
        #classname = make_classname(prefix, "finite_element", ir["element_numbers"][ufl_element])


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

        table_types = self.ir["uflacs"]["expr_irs"][num_points]["table_types"]
        tt = table_types[uname]
        if tt == "zeros":
            error("Not expecting zero arguments to get this far.")
            return L.LiteralFloat(0.0)
        elif tt == "ones":
            warning("Should simplify ones arguments before getting this far.")
            return L.LiteralFloat(1.0)

        entity = self.symbols.entity(self.ir["entitytype"], mt.restriction)
        idof = self.symbols.argument_loop_index(mt.terminal.number())

        if tt == "piecewise":
            iq = 0
        else:
            iq = self.symbols.quadrature_loop_index(num_points)

        uname = L.Symbol(uname)
        return uname[entity][iq][idof - begin]


    def coefficient(self, e, mt, tabledata, num_points):
        # TODO: Passing type along with tabledata would make a lot of code cleaner
        uname, begin, end = tabledata
        table_types = self.ir["uflacs"]["expr_irs"][num_points]["table_types"]
        tt = table_types[uname]

        if tt == "zeros":
            # FIXME: Remove at earlier stage so dependent code can also be removed
            warning("Not expecting zero coefficients to get this far.")
            L = self.language
            return L.LiteralFloat(0.0)
        elif tt == "ones" and (end - begin) == 1:
            # f = 1.0 * f_i, just return direct reference to dof array at dof begin
            return self.symbols.coefficient_dof_access(mt.terminal, begin)
        else:
            # Return symbol, see definitions for computation 
            return self.symbols.coefficient_value(mt)  #, num_points)


    def quadrature_weight(self, e, mt, tabledata, num_points):
        weight = self.symbols.weights_array(num_points)
        iq = self.symbols.quadrature_loop_index(num_points)
        return weight[iq]


    def spatial_coordinate(self, e, mt, tabledata, num_points):
        #L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of SpatialCoordinate.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of SpatialCoordinate.")
        ffc_assert(not mt.averaged, "Not expecting average of SpatialCoordinates.")

        if self.ir["integral_type"] in physical_quadrature_integral_types:
            # Physical coordinates are available in given variables
            assert num_points is None
            x = self.symbols.points_array(num_points)
            iq = self.symbols.quadrature_loop_index(num_points)
            gdim, = mt.terminal.ufl_shape
            return x[iq * gdim + mt.flat_component]
        else:
            # Physical coordinates are computed by code generated in definitions
            return self.symbols.x_component(mt)


    def cell_coordinate(self, e, mt, tabledata, num_points):
        #L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of CellCoordinate.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of CellCoordinate.")
        ffc_assert(not mt.averaged, "Not expecting average of CellCoordinate.")

        if self.ir["integral_type"] == "cell" and not mt.restriction:
            X = self.symbols.points_array(num_points)
            if num_points == 1:
                return X[mt.flat_component]
            else:
                iq = self.symbols.quadrature_loop_index(num_points)
                tdim, = mt.terminal.ufl_shape
                return X[iq * tdim + mt.flat_component]
        else:
            # X should be computed from x or Xf symbolically instead of getting here
            error("Expecting reference cell coordinate to be symbolically rewritten.")


    def facet_coordinate(self, e, mt, tabledata, num_points):
        L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of FacetCoordinate.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of FacetCoordinate.")
        ffc_assert(not mt.averaged, "Not expecting average of FacetCoordinate.")
        ffc_assert(not mt.restriction, "Not expecting restriction of FacetCoordinate.")

        if self.ir["integral_type"] in ("interior_facet", "exterior_facet"):
            tdim, = mt.terminal.ufl_shape
            if tdim == 0:
                error("Vertices have no facet coordinates.")
            elif tdim == 1:
                # 0D vertex coordinate
                warning("Vertex coordinate is always 0, should get rid of this in ufl geometry lowering.")
                return L.LiteralFloat(0.0)

            Xf = self.points_array(num_points)
            iq = self.symbols.quadrature_loop_index(num_points)
            assert 0 <= mt.flat_component < (tdim-1)
            if tdim == 2:
                # 1D edge coordinate
                assert mt.flat_component == 0
                return Xf[iq]
            else:
                # The general case
                return Xf[iq * (tdim - 1) + mt.flat_component]
        else:
            # Xf should be computed from X or x symbolically instead of getting here
            error("Expecting reference facet coordinate to be symbolically rewritten.")


    def jacobian(self, e, mt, tabledata, num_points):
        L = self.language
        ffc_assert(not mt.global_derivatives, "Not expecting derivatives of Jacobian.")
        ffc_assert(not mt.local_derivatives, "Not expecting derivatives of Jacobian.")
        ffc_assert(not mt.averaged, "Not expecting average of Jacobian.")
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
