# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFC/UFC specific variable access."""

import logging
import warnings

import ufl
from ffc.fiatinterface import create_element
from ufl.finiteelement import MixedElement
from ufl.measure import custom_integral_types

logger = logging.getLogger(__name__)


class FFCBackendAccess(object):
    """FFC specific cpp formatter class."""

    def __init__(self, ir, language, symbols, parameters):

        # Store ir and parameters
        self.entitytype = ir.entitytype
        self.integral_type = ir.integral_type
        self.language = language
        self.symbols = symbols
        self.parameters = parameters

        # Lookup table for handler to call when the "get" method (below) is
        # called, depending on the first argument type.
        self.call_lookup = {ufl.coefficient.Coefficient: self.coefficient,
                            ufl.geometry.Jacobian: self.jacobian,
                            ufl.geometry.CellCoordinate: self.cell_coordinate,
                            ufl.geometry.FacetCoordinate: self.facet_coordinate,
                            ufl.geometry.CellVertices: self.cell_vertices,
                            ufl.geometry.FacetEdgeVectors: self.facet_edge_vectors,
                            ufl.geometry.CellEdgeVectors: self.cell_edge_vectors,
                            ufl.geometry.CellFacetJacobian: self.cell_facet_jacobian,
                            ufl.geometry.ReferenceCellVolume: self.reference_cell_volume,
                            ufl.geometry.ReferenceFacetVolume: self.reference_facet_volume,
                            ufl.geometry.ReferenceCellEdgeVectors: self.reference_cell_edge_vectors,
                            ufl.geometry.ReferenceFacetEdgeVectors: self.reference_facet_edge_vectors,
                            ufl.geometry.ReferenceNormal: self.reference_normal,
                            ufl.geometry.FacetOrientation: self.facet_orientation,
                            ufl.geometry.CellOrientation: self.cell_orientation,
                            ufl.geometry.SpatialCoordinate: self.spatial_coordinate}

    def get(self, e, mt, tabledata, num_points):
        # Call appropriate handler, depending on the type of e
        etype = type(e)
        if etype in self.call_lookup:
            return self.call_lookup[etype](e, mt, tabledata, num_points)
        else:
            raise RuntimeError("Not handled: %s", etype)

    def coefficient(self, e, mt, tabledata, num_points):
        ttype = tabledata.ttype
        assert ttype != "zeros"

        begin, end = tabledata.dofrange

        if ttype == "ones" and (end - begin) == 1:
            # f = 1.0 * f_{begin}, just return direct reference to dof
            # array at dof begin (if mt is restricted, begin contains
            # cell offset)
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
            return self.symbols.coefficient_value(mt)  # , num_points)

    def spatial_coordinate(self, e, mt, tabledata, num_points):
        if mt.global_derivatives:
            raise RuntimeError("Not expecting global derivatives of SpatialCoordinate.")
        if mt.averaged:
            raise RuntimeError("Not expecting average of SpatialCoordinates.")

        if self.integral_type in custom_integral_types:
            if mt.local_derivatives:
                raise RuntimeError("FIXME: Jacobian in custom integrals is not implemented.")

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
            # Physical coordinates are computed by code generated in
            # definitions
            return self.symbols.x_component(mt)

    def cell_coordinate(self, e, mt, tabledata, num_points):
        if mt.global_derivatives:
            raise RuntimeError("Not expecting derivatives of CellCoordinate.")
        if mt.local_derivatives:
            raise RuntimeError("Not expecting derivatives of CellCoordinate.")
        if mt.averaged:
            raise RuntimeError("Not expecting average of CellCoordinate.")

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
            # X should be computed from x or Xf symbolically instead of
            # getting here
            raise RuntimeError("Expecting reference cell coordinate to be symbolically rewritten.")

    def facet_coordinate(self, e, mt, tabledata, num_points):
        L = self.language
        if mt.global_derivatives:
            raise RuntimeError("Not expecting derivatives of FacetCoordinate.")
        if mt.local_derivatives:
            raise RuntimeError("Not expecting derivatives of FacetCoordinate.")
        if mt.averaged:
            raise RuntimeError("Not expecting average of FacetCoordinate.")
        if mt.restriction:
            raise RuntimeError("Not expecting restriction of FacetCoordinate.")

        if self.integral_type in ("interior_facet", "exterior_facet"):
            tdim, = mt.terminal.ufl_shape
            if tdim == 0:
                raise RuntimeError("Vertices have no facet coordinates.")
            elif tdim == 1:
                warnings.warn(
                    "Vertex coordinate is always 0, should get rid of this in ufl geometry lowering."
                )
                return L.LiteralFloat(0.0)
            Xf = self.points_table(num_points)
            iq = self.symbols.quadrature_loop_index()
            assert 0 <= mt.flat_component < (tdim - 1)
            if num_points == 1:
                index = mt.flat_component
            elif tdim == 2:
                index = iq
            else:
                index = iq * (tdim - 1) + mt.flat_component
            return Xf[index]
        else:
            # Xf should be computed from X or x symbolically instead of
            # getting here
            raise RuntimeError("Expecting reference facet coordinate to be symbolically rewritten.")

    def jacobian(self, e, mt, tabledata, num_points):
        if mt.global_derivatives:
            raise RuntimeError("Not expecting global derivatives of Jacobian.")
        if mt.averaged:
            raise RuntimeError("Not expecting average of Jacobian.")
        return self.symbols.J_component(mt)

    def reference_cell_volume(self, e, mt, tabledata, access):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            return L.Symbol("{0}_reference_cell_volume".format(cellname))
        else:
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

    def reference_facet_volume(self, e, mt, tabledata, access):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            return L.Symbol("{0}_reference_facet_volume".format(cellname))
        else:
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

    def reference_normal(self, e, mt, tabledata, access):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol("{0}_reference_facet_normals".format(cellname))
            facet = self.symbols.entity("facet", mt.restriction)
            return table[facet][mt.component[0]]
        else:
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

    def cell_facet_jacobian(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol("{0}_reference_facet_jacobian".format(cellname))
            facet = self.symbols.entity("facet", mt.restriction)
            return table[facet][mt.component[0]][mt.component[1]]
        elif cellname == "interval":
            raise RuntimeError("The reference facet jacobian doesn't make sense for interval cell.")
        else:
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

    def reference_cell_edge_vectors(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol("{0}_reference_edge_vectors".format(cellname))
            return table[mt.component[0]][mt.component[1]]
        elif cellname == "interval":
            raise RuntimeError("The reference cell edge vectors doesn't make sense for interval cell.")
        else:
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

    def reference_facet_edge_vectors(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname in ("tetrahedron", "hexahedron"):
            table = L.Symbol("{0}_reference_edge_vectors".format(cellname))
            facet = self.symbols.entity("facet", mt.restriction)
            return table[facet][mt.component[0]][mt.component[1]]
        elif cellname in ("interval", "triangle", "quadrilateral"):
            raise RuntimeError(
                "The reference cell facet edge vectors doesn't make sense for interval or triangle cell."
            )
        else:
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

    def cell_orientation(self, e, mt, tabledata, num_points):
        # Error if not in manifold case:
        domain = mt.terminal.ufl_domain()
        assert domain.geometric_dimension() > domain.topological_dimension()
        return self.symbols.cell_orientation_internal(mt.restriction)

    def facet_orientation(self, e, mt, tabledata, num_points):
        L = self.language
        cellname = mt.terminal.ufl_domain().ufl_cell().cellname()
        if cellname not in ("interval", "triangle", "tetrahedron"):
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

        table = L.Symbol("{0}_facet_orientations".format(cellname))
        facet = self.symbols.entity("facet", mt.restriction)
        return table[facet]

    def cell_vertices(self, e, mt, tabledata, num_points):
        # Get properties of domain
        domain = mt.terminal.ufl_domain()
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()

        # Get dimension and dofmap of scalar element
        assert isinstance(coordinate_element, MixedElement)
        assert coordinate_element.value_shape() == (gdim, )
        ufl_scalar_element, = set(coordinate_element.sub_elements())
        assert ufl_scalar_element.family() in ("Lagrange", "Q", "S")
        fiat_scalar_element = create_element(ufl_scalar_element)
        vertex_scalar_dofs = fiat_scalar_element.entity_dofs()[0]
        num_scalar_dofs = fiat_scalar_element.space_dimension()

        # Get dof and component
        dof, = vertex_scalar_dofs[mt.component[0]]
        component = mt.component[1]

        expr = self.symbols.domain_dof_access(dof, component, gdim, num_scalar_dofs, mt.restriction)
        return expr

    def cell_edge_vectors(self, e, mt, tabledata, num_points):
        # Get properties of domain
        domain = mt.terminal.ufl_domain()
        cellname = domain.ufl_cell().cellname()
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()

        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            pass
        elif cellname == "interval":
            raise RuntimeError("The physical cell edge vectors doesn't make sense for interval cell.")
        else:
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

        # Get dimension and dofmap of scalar element
        assert isinstance(coordinate_element, MixedElement)
        assert coordinate_element.value_shape() == (gdim, )
        ufl_scalar_element, = set(coordinate_element.sub_elements())
        assert ufl_scalar_element.family() in ("Lagrange", "Q", "S")
        fiat_scalar_element = create_element(ufl_scalar_element)
        vertex_scalar_dofs = fiat_scalar_element.entity_dofs()[0]
        num_scalar_dofs = fiat_scalar_element.space_dimension()

        # Get edge vertices
        edge = mt.component[0]
        edge_vertices = fiat_scalar_element.get_reference_element().get_topology()[1][edge]
        vertex0, vertex1 = edge_vertices

        # Get dofs and component
        dof0, = vertex_scalar_dofs[vertex0]
        dof1, = vertex_scalar_dofs[vertex1]
        component = mt.component[1]

        expr = (self.symbols.domain_dof_access(dof0, component,
                                               gdim, num_scalar_dofs,
                                               mt.restriction)
                - self.symbols.domain_dof_access(dof1, component,
                                                 gdim, num_scalar_dofs,
                                                 mt.restriction))
        return expr

    def facet_edge_vectors(self, e, mt, tabledata, num_points):
        L = self.language

        # Get properties of domain
        domain = mt.terminal.ufl_domain()
        cellname = domain.ufl_cell().cellname()
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()

        if cellname in ("tetrahedron", "hexahedron"):
            pass
        elif cellname in ("interval", "triangle", "quadrilateral"):
            raise RuntimeError(
                "The physical facet edge vectors doesn't make sense for {0} cell.".format(cellname))
        else:
            raise RuntimeError("Unhandled cell types {0}.".format(cellname))

        # Get dimension and dofmap of scalar element
        assert isinstance(coordinate_element, MixedElement)
        assert coordinate_element.value_shape() == (gdim, )
        ufl_scalar_element, = set(coordinate_element.sub_elements())
        assert ufl_scalar_element.family() in ("Lagrange", "Q", "S")
        fiat_scalar_element = create_element(ufl_scalar_element)
        num_scalar_dofs = fiat_scalar_element.space_dimension()

        # Get edge vertices
        facet = self.symbols.entity("facet", mt.restriction)
        facet_edge = mt.component[0]
        facet_edge_vertices = L.Symbol("{0}_facet_edge_vertices".format(cellname))
        vertex0 = facet_edge_vertices[facet][facet_edge][0]
        vertex1 = facet_edge_vertices[facet][facet_edge][1]

        # Get dofs and component
        component = mt.component[1]
        assert coordinate_element.degree() == 1, "Assuming degree 1 element"
        dof0 = vertex0
        dof1 = vertex1
        expr = (
            self.symbols.domain_dof_access(dof0, component, gdim, num_scalar_dofs, mt.restriction)
            - self.symbols.domain_dof_access(dof1, component, gdim, num_scalar_dofs, mt.restriction))

        return expr
