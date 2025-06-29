# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCx/UFC specific variable access."""

import logging
import warnings
from typing import Optional

import basix.ufl
import ufl

import ffcx.codegeneration.lnodes as L
from ffcx.definitions import entity_types
from ffcx.ir.analysis.modified_terminals import ModifiedTerminal
from ffcx.ir.elementtables import UniqueTableReferenceT
from ffcx.ir.representationutils import QuadratureRule

logger = logging.getLogger("ffcx")


class FFCXBackendAccess:
    """FFCx specific formatter class."""

    entity_type: entity_types

    def __init__(self, entity_type: entity_types, integral_type: str, symbols, options):
        """Initialise."""
        # Store ir and options
        self.entity_type = entity_type
        self.integral_type = integral_type
        self.symbols = symbols
        self.options = options

        # Lookup table for handler to call when the "get" method (below) is
        # called, depending on the first argument type.
        self.call_lookup = {
            ufl.coefficient.Coefficient: self.coefficient,
            ufl.constant.Constant: self.constant,
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
            ufl.geometry.CellOrientation: self._pass,
            ufl.geometry.FacetOrientation: self.facet_orientation,
            ufl.geometry.SpatialCoordinate: self.spatial_coordinate,
        }

    def get(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
    ):
        """Format a terminal."""
        e = mt.terminal
        # Call appropriate handler, depending on the type of e
        handler = self.call_lookup.get(type(e), False)  # type: ignore

        if not handler:
            # Look for parent class types instead
            for k in self.call_lookup.keys():
                if isinstance(e, k):
                    handler = self.call_lookup[k]
                    break

        if handler:
            return handler(mt, tabledata, quadrature_rule)
        else:
            raise RuntimeError(f"Not handled: {type(e)}")

    def coefficient(
        self,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule,
    ):
        """Access a coefficient."""
        ttype = tabledata.ttype
        assert ttype != "zeros"

        num_dofs = tabledata.values.shape[3]
        begin = tabledata.offset
        assert begin is not None
        assert tabledata.block_size is not None
        end = begin + tabledata.block_size * (num_dofs - 1) + 1

        if ttype == "ones" and (end - begin) == 1:
            # f = 1.0 * f_{begin}, just return direct reference to dof
            # array at dof begin (if mt is restricted, begin contains
            # cell offset)
            return self.symbols.coefficient_dof_access(mt.terminal, begin)
        else:
            # Return symbol, see definitions for computation
            return self.symbols.coefficient_value(mt)

    def constant(
        self,
        mt: ModifiedTerminal,
        tabledata: Optional[UniqueTableReferenceT],
        quadrature_rule: Optional[QuadratureRule],
    ):
        """Access a constant."""
        # Access to a constant is handled trivially, directly through constants symbol
        return self.symbols.constant_index_access(mt.terminal, mt.flat_component)

    def spatial_coordinate(
        self, mt: ModifiedTerminal, tabledata: UniqueTableReferenceT, num_points: QuadratureRule
    ):
        """Access a spatial coordinate."""
        if mt.global_derivatives:
            raise RuntimeError("Not expecting global derivatives of SpatialCoordinate.")
        if mt.averaged is not None:
            raise RuntimeError("Not expecting average of SpatialCoordinates.")

        if self.integral_type in ufl.custom_integral_types:
            if mt.local_derivatives:
                raise RuntimeError("FIXME: Jacobian in custom integrals is not implemented.")

            # Access predefined quadrature points table
            x = self.symbols.custom_points_table
            iq = self.symbols.quadrature_loop_index
            (gdim,) = mt.terminal.ufl_shape
            if gdim == 1:
                index = iq
            else:
                index = iq * gdim + mt.flat_component
            return x[index]
        elif self.integral_type == "expression":
            # Physical coordinates are computed by code generated in
            # definitions
            return self.symbols.x_component(mt)
        else:
            # Physical coordinates are computed by code generated in
            # definitions
            return self.symbols.x_component(mt)

    def cell_coordinate(self, mt, tabledata, num_points):
        """Access a cell coordinate."""
        if mt.global_derivatives:
            raise RuntimeError("Not expecting derivatives of CellCoordinate.")
        if mt.local_derivatives:
            raise RuntimeError("Not expecting derivatives of CellCoordinate.")
        if mt.averaged is not None:
            raise RuntimeError("Not expecting average of CellCoordinate.")

        if self.integral_type == "cell" and not mt.restriction:
            # Access predefined quadrature points table
            X = self.symbols.points_table(num_points)
            (tdim,) = mt.terminal.ufl_shape
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

    def facet_coordinate(self, mt, tabledata, num_points):
        """Access a facet coordinate."""
        if mt.global_derivatives:
            raise RuntimeError("Not expecting derivatives of FacetCoordinate.")
        if mt.local_derivatives:
            raise RuntimeError("Not expecting derivatives of FacetCoordinate.")
        if mt.averaged is not None:
            raise RuntimeError("Not expecting average of FacetCoordinate.")
        if mt.restriction:
            raise RuntimeError("Not expecting restriction of FacetCoordinate.")

        if self.integral_type in ("interior_facet", "exterior_facet"):
            (tdim,) = mt.terminal.ufl_shape
            if tdim == 0:
                raise RuntimeError("Vertices have no facet coordinates.")
            elif tdim == 1:
                warnings.warn(
                    "Vertex coordinate is always 0, should get rid of this in UFL "
                    "geometry lowering."
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

    def jacobian(self, mt, tabledata, num_points):
        """Access a jacobian."""
        if mt.averaged is not None:
            raise RuntimeError("Not expecting average of Jacobian.")
        return self.symbols.J_component(mt)

    def reference_cell_volume(self, mt, tabledata, access):
        """Access a reference cell volume."""
        cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            return L.Symbol(f"{cellname}_reference_cell_volume", dtype=L.DataType.REAL)
        else:
            raise RuntimeError(f"Unhandled cell types {cellname}.")

    def reference_facet_volume(self, mt, tabledata, access):
        """Access a reference facet volume."""
        cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            return L.Symbol(f"{cellname}_reference_facet_volume", dtype=L.DataType.REAL)
        else:
            raise RuntimeError(f"Unhandled cell types {cellname}.")

    def reference_normal(self, mt, tabledata, access):
        """Access a reference normal."""
        cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
        if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol(f"{cellname}_reference_normals", dtype=L.DataType.REAL)
            facet = self.symbols.entity("facet", mt.restriction)
            return table[facet][mt.component[0]]
        else:
            raise RuntimeError(f"Unhandled cell types {cellname}.")

    def cell_facet_jacobian(self, mt, tabledata, num_points):
        """Access a cell facet jacobian."""
        cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
        if cellname in (
            "triangle",
            "tetrahedron",
            "quadrilateral",
            "hexahedron",
            "prism",
            "pyramid",
        ):
            table = L.Symbol(f"{cellname}_cell_facet_jacobian", dtype=L.DataType.REAL)
            facet = self.symbols.entity("facet", mt.restriction)
            return table[facet][mt.component[0]][mt.component[1]]
        elif cellname == "interval":
            raise RuntimeError("The reference facet jacobian doesn't make sense for interval cell.")
        else:
            raise RuntimeError(f"Unhandled cell types {cellname}.")

    def reference_cell_edge_vectors(self, mt, tabledata, num_points):
        """Access a reference cell edge vector."""
        cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            table = L.Symbol(f"{cellname}_reference_cell_edge_vectors", dtype=L.DataType.REAL)
            return table[mt.component[0]][mt.component[1]]
        elif cellname == "interval":
            raise RuntimeError(
                "The reference cell edge vectors doesn't make sense for interval cell."
            )
        else:
            raise RuntimeError(f"Unhandled cell types {cellname}.")

    def reference_facet_edge_vectors(self, mt, tabledata, num_points):
        """Access a reference facet edge vector."""
        cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
        if cellname in ("tetrahedron", "hexahedron"):
            table = L.Symbol(f"{cellname}_reference_facet_edge_vectors", dtype=L.DataType.REAL)
            return table[mt.component[0]][mt.component[1]]
        elif cellname in ("interval", "triangle", "quadrilateral"):
            raise RuntimeError(
                "The reference cell facet edge vectors doesn't make sense for interval "
                "or triangle cell."
            )
        else:
            raise RuntimeError(f"Unhandled cell types {cellname}.")

    def facet_orientation(self, mt, tabledata, num_points):
        """Access a facet orientation."""
        cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
        if cellname not in ("interval", "triangle", "tetrahedron"):
            raise RuntimeError(f"Unhandled cell types {cellname}.")

        table = L.Symbol(f"{cellname}_facet_orientation", dtype=L.DataType.INT)
        facet = self.symbols.entity("facet", mt.restriction)
        return table[facet]

    def cell_vertices(self, mt, tabledata, num_points):
        """Access a cell vertex."""
        # Get properties of domain
        domain = ufl.domain.extract_unique_domain(mt.terminal)
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()

        # Get dimension and dofmap of scalar element
        assert isinstance(coordinate_element, basix.ufl._BlockedElement)
        assert coordinate_element.reference_value_shape == (gdim,)
        (ufl_scalar_element,) = set(coordinate_element.sub_elements)
        scalar_element = ufl_scalar_element
        assert scalar_element.reference_value_size == 1 and scalar_element.block_size == 1

        vertex_scalar_dofs = scalar_element.entity_dofs[0]
        num_scalar_dofs = scalar_element.dim

        # Get dof and component
        (dof,) = vertex_scalar_dofs[mt.component[0]]
        component = mt.component[1]

        expr = self.symbols.domain_dof_access(dof, component, gdim, num_scalar_dofs, mt.restriction)
        return expr

    def cell_edge_vectors(self, mt, tabledata, num_points):
        """Access a cell edge vector."""
        # Get properties of domain
        domain = ufl.domain.extract_unique_domain(mt.terminal)
        cellname = domain.ufl_cell().cellname()
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()

        if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
            pass
        elif cellname == "interval":
            raise RuntimeError(
                "The physical cell edge vectors doesn't make sense for interval cell."
            )
        else:
            raise RuntimeError(f"Unhandled cell types {cellname}.")

        # Get dimension and dofmap of scalar element
        assert isinstance(coordinate_element, basix.ufl._BlockedElement)
        assert coordinate_element.reference_value_shape == (gdim,)
        (ufl_scalar_element,) = set(coordinate_element.sub_elements)
        scalar_element = ufl_scalar_element
        assert scalar_element.reference_value_size == 1 and scalar_element.block_size == 1

        vertex_scalar_dofs = scalar_element.entity_dofs[0]
        num_scalar_dofs = scalar_element.dim

        # Get edge vertices
        edge = mt.component[0]
        vertex0, vertex1 = scalar_element.reference_topology[1][edge]

        # Get dofs and component
        (dof0,) = vertex_scalar_dofs[vertex0]
        (dof1,) = vertex_scalar_dofs[vertex1]
        component = mt.component[1]

        return self.symbols.domain_dof_access(
            dof0, component, gdim, num_scalar_dofs, mt.restriction
        ) - self.symbols.domain_dof_access(dof1, component, gdim, num_scalar_dofs, mt.restriction)

    def facet_edge_vectors(self, mt, tabledata, num_points):
        """Access a facet edge vector."""
        # Get properties of domain
        domain = ufl.domain.extract_unique_domain(mt.terminal)
        cellname = domain.ufl_cell().cellname()
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()

        if cellname in ("tetrahedron", "hexahedron"):
            pass
        elif cellname in ("interval", "triangle", "quadrilateral"):
            raise RuntimeError(
                f"The physical facet edge vectors doesn't make sense for {cellname} cell."
            )
        else:
            raise RuntimeError(f"Unhandled cell types {cellname}.")

        # Get dimension and dofmap of scalar element
        assert isinstance(coordinate_element, basix.ufl._BlockedElement)
        assert coordinate_element.reference_value_shape == (gdim,)
        (ufl_scalar_element,) = set(coordinate_element.sub_elements)
        scalar_element = ufl_scalar_element
        assert scalar_element.reference_value_size == 1 and scalar_element.block_size == 1

        scalar_element = ufl_scalar_element
        num_scalar_dofs = scalar_element.dim

        # Get edge vertices
        facet = self.symbols.entity("facet", mt.restriction)
        facet_edge = mt.component[0]
        facet_edge_vertices = L.Symbol(f"{cellname}_facet_edge_vertices", dtype=L.DataType.INT)
        vertex0 = facet_edge_vertices[facet][facet_edge][0]
        vertex1 = facet_edge_vertices[facet][facet_edge][1]

        # Get dofs and component
        component = mt.component[1]
        assert coordinate_element.embedded_superdegree == 1, "Assuming degree 1 element"
        dof0 = vertex0
        dof1 = vertex1
        expr = self.symbols.domain_dof_access(
            dof0, component, gdim, num_scalar_dofs, mt.restriction
        ) - self.symbols.domain_dof_access(dof1, component, gdim, num_scalar_dofs, mt.restriction)

        return expr

    def _pass(self, *args, **kwargs):
        """Return one."""
        return 1

    def table_access(
        self,
        tabledata: UniqueTableReferenceT,
        entity_type: entity_types,
        restriction: str,
        quadrature_index: L.MultiIndex,
        dof_index: L.MultiIndex,
    ):
        """Access element table for given entity, quadrature point, and dof index.

        Args:
            tabledata: Table data object
            entity_type: Entity type
            restriction: Restriction ("+", "-")
            quadrature_index: Quadrature index
            dof_index: Dof index
        """
        entity = self.symbols.entity(entity_type, restriction)
        iq_global_index = quadrature_index.global_index
        ic_global_index = dof_index.global_index
        qp = 0  # quadrature permutation

        symbols = []
        if tabledata.is_uniform:
            entity = L.LiteralInt(0)

        if tabledata.is_piecewise:
            iq_global_index = L.LiteralInt(0)

        # FIXME: Hopefully tabledata is not permuted when applying sum
        # factorization
        if tabledata.is_permuted:
            qp = self.symbols.quadrature_permutation[0]
            if restriction == "-":
                qp = self.symbols.quadrature_permutation[1]

        if dof_index.dim == 1 and quadrature_index.dim == 1:
            symbols += [L.Symbol(tabledata.name, dtype=L.DataType.REAL)]
            return self.symbols.element_tables[tabledata.name][qp][entity][iq_global_index][
                ic_global_index
            ], symbols
        else:
            FE = []
            assert tabledata.tensor_factors is not None
            for i in range(dof_index.dim):
                factor = tabledata.tensor_factors[i]
                iq_i = quadrature_index.local_index(i)
                ic_i = dof_index.local_index(i)
                table = self.symbols.element_tables[factor.name][qp][entity][iq_i][ic_i]
                symbols += [L.Symbol(factor.name, dtype=L.DataType.REAL)]
                FE.append(table)
            return L.Product(FE), symbols
