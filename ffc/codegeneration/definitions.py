# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFC/UFC specific variable definitions."""

import logging

import ufl

logger = logging.getLogger(__name__)


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


class FFCBackendDefinitions(object):
    """FFC specific code definitions."""

    def __init__(self, ir, language, symbols, parameters):
        # Store ir and parameters
        self.integral_type = ir.integral_type
        self.entitytype = ir.entitytype
        self.language = language
        self.symbols = symbols
        self.parameters = parameters

        # Lookup table for handler to call when the "get" method (below) is
        # called, depending on the first argument type.
        self.call_lookup = {ufl.coefficient.Coefficient: self.coefficient,
                            ufl.geometry.Jacobian: self.jacobian,
                            ufl.geometry.CellVertices: self._expect_physical_coords,
                            ufl.geometry.FacetEdgeVectors: self._expect_physical_coords,
                            ufl.geometry.CellEdgeVectors: self._expect_physical_coords,
                            ufl.geometry.CellFacetJacobian: self._expect_table,
                            ufl.geometry.ReferenceCellVolume: self._expect_table,
                            ufl.geometry.ReferenceFacetVolume: self._expect_table,
                            ufl.geometry.ReferenceCellEdgeVectors: self._expect_table,
                            ufl.geometry.ReferenceFacetEdgeVectors: self._expect_table,
                            ufl.geometry.ReferenceNormal: self._expect_table,
                            ufl.geometry.FacetOrientation: self._expect_table,
                            ufl.geometry.CellOrientation: self.cell_orientation,
                            ufl.geometry.SpatialCoordinate: self.spatial_coordinate}

    def get(self, t, mt, tabledata, num_points, access):
        # Call appropriate handler, depending on the type of t
        ttype = type(t)
        if ttype in self.call_lookup:
            return self.call_lookup[ttype](t, mt, tabledata, num_points, access)
        else:
            raise RuntimeError("Not handled: %s", ttype)

    def coefficient(self, t, mt, tabledata, num_points, access):
        """Return definition code for coefficients."""
        L = self.language

        ttype = tabledata.ttype
        begin, end = tabledata.dofrange

        # fe_classname = ir.classnames["finite_element"][t.ufl_element()]

        if ttype == "zeros":
            logging.debug("Not expecting zero coefficients to get this far.")
            return []

        # For a constant coefficient we reference the dofs directly, so no definition needed
        if ttype == "ones" and (end - begin) == 1:
            return []

        # For quadrature elements we reference the dofs directly, so no definition needed
        if ttype == "quadrature":
            return []

        assert begin < end

        # Get access to element table
        FE = self.symbols.element_table(tabledata, self.entitytype, mt.restriction)

        unroll = len(tabledata.dofmap) != end - begin
        # unroll = True
        if unroll:
            # TODO: Could also use a generated constant dofmap here like in block code
            # Unrolled loop to accumulate linear combination of dofs and tables
            values = [
                self.symbols.coefficient_dof_access(mt.terminal, idof) * FE[i]
                for i, idof in enumerate(tabledata.dofmap)
            ]
            value = L.Sum(values)
            code = [L.VariableDecl("const ufc_scalar_t", access, value)]
        else:
            # Loop to accumulate linear combination of dofs and tables
            ic = self.symbols.coefficient_dof_sum_index()
            dof_access = self.symbols.coefficient_dof_access(mt.terminal, ic + begin)
            code = [
                L.VariableDecl("ufc_scalar_t", access, 0.0),
                L.ForRange(ic, 0, end - begin, body=[L.AssignAdd(access, dof_access * FE[ic])])
            ]
        return code

    def _define_coordinate_dofs_lincomb(self, e, mt, tabledata, num_points, access):
        """Define x or J as a linear combination of coordinate dofs with given table data."""
        L = self.language

        # Get properties of domain
        domain = mt.terminal.ufl_domain()
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()
        num_scalar_dofs = num_coordinate_component_dofs(coordinate_element)

        # Reference coordinates are known, no coordinate field, so we compute
        # this component as linear combination of coordinate_dofs "dofs" and table

        # Find table name and dof range it corresponds to
        ttype = tabledata.ttype
        begin, end = tabledata.dofrange

        assert end - begin <= num_scalar_dofs
        assert ttype != "zeros"
        assert ttype != "ones"
        assert ttype != "quadrature"

        # Get access to element table
        #        if tabledata.is_uniform and tabledata.is_piecewise:
        #            FE = tabledata.values[0][0]
        #        else:
        FE = self.symbols.element_table(tabledata, self.entitytype, mt.restriction)

        # Inlined version (we know this is bounded by a small number)
        dof_access = self.symbols.domain_dofs_access(gdim, num_scalar_dofs, mt.restriction)
        value = L.Sum([dof_access[idof] * FE[i] for i, idof in enumerate(tabledata.dofmap)])
        code = [L.VariableDecl("const double", access, value)]

        return code

    def spatial_coordinate(self, e, mt, tabledata, num_points, access):
        """Return definition code for the physical spatial coordinates.

        If physical coordinates are given:
          No definition needed.

        If reference coordinates are given:
          x = sum_k xdof_k xphi_k(X)

        If reference facet coordinates are given:
          x = sum_k xdof_k xphi_k(Xf)
        """
        if self.integral_type in ufl.measure.custom_integral_types:
            # FIXME: Jacobian may need adjustment for custom_integral_types
            if mt.local_derivatives:
                logging.exception("FIXME: Jacobian in custom integrals is not implemented.")
            return []
        else:
            return self._define_coordinate_dofs_lincomb(e, mt, tabledata, num_points, access)

    def jacobian(self, e, mt, tabledata, num_points, access):
        """Return definition code for the Jacobian of x(X).

        J = sum_k xdof_k grad_X xphi_k(X)
        """
        # TODO: Jacobian may need adjustment for custom_integral_types
        return self._define_coordinate_dofs_lincomb(e, mt, tabledata, num_points, access)

    def cell_orientation(self, e, mt, tabledata, num_points, access):
        # Would be nicer if cell_orientation was a double variable input,
        # but this is how dolfin/ufc/ffc currently passes this information.
        # 0 means up and gives +1.0, 1 means down and gives -1.0.
        L = self.language
        co = self.symbols.cell_orientation_argument(mt.restriction)
        expr = L.Conditional(L.EQ(co, L.LiteralInt(1)), L.LiteralFloat(-1.0), L.LiteralFloat(+1.0))
        code = [L.VariableDecl("const ufc_scalar_t", access, expr)]
        return code

    def _expect_table(self, e, mt, tabledata, num_points, access):
        """These quantities refer to constant tables defined in ufc_geometry.h."""
        # TODO: Inject const static table here instead?
        return []

    def _expect_physical_coords(self, e, mt, tabledata, num_points, access):
        """These quantities refer to coordinate_dofs"""
        # TODO: Generate more efficient inline code for Max/MinCell/FacetEdgeLength
        #       and CellDiameter here rather than lowering these quantities?
        return []
