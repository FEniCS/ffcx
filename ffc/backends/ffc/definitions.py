# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFC/UFC specific variable definitions."""

import logging

import ufl
from ffc.backends.ffc.common import num_coordinate_component_dofs

logger = logging.getLogger(__name__)


class FFCBackendDefinitions(ufl.corealg.multifunction.MultiFunction):
    """FFC specific code definitions."""

    def __init__(self, ir, language, symbols, parameters):
        ufl.corealg.multifunction.MultiFunction.__init__(self)

        # Store ir and parameters
        self.integral_type = ir["integral_type"]
        self.entitytype = ir["entitytype"]
        self.language = language
        self.symbols = symbols
        self.parameters = parameters

    # === Generate code to define variables for ufl types ===

    def expr(self, t, mt, tabledata, num_points, access):
        logging.exception("Unhandled type {0}".format(type(t)))

    # def quadrature_weight(self, e, mt, tabledata, num_points, access):
    #    "Quadrature weights are precomputed and need no code."
    #    return []

    def constant_value(self, e, mt, tabledata, num_points, access):
        """Constants simply use literals in the target language."""
        return []

    def argument(self, t, mt, tabledata, num_points, access):
        """Arguments are accessed through element tables."""
        return []

    def coefficient(self, t, mt, tabledata, num_points, access):
        """Return definition code for coefficients."""
        L = self.language

        ttype = tabledata.ttype
        begin, end = tabledata.dofrange

        # fe_classname = ir["classnames"]["finite_element"][t.ufl_element()]

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
            code = [L.VariableDecl("const double", access, value)]
        else:
            # Loop to accumulate linear combination of dofs and tables
            ic = self.symbols.coefficient_dof_sum_index()
            dof_access = self.symbols.coefficient_dof_access(mt.terminal, ic + begin)
            code = [
                L.VariableDecl("double", access, 0.0),
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

    def cell_coordinate(self, e, mt, tabledata, num_points, access):
        """Return definition code for the reference spatial coordinates.

        If reference coordinates are given::

            No definition needed.

        If physical coordinates are given and domain is affine::

            X = K*(x-x0)

        This is inserted symbolically.

        If physical coordinates are given and domain is non- affine::

            Not currently supported.

        """
        # Should be either direct access to points array or symbolically computed
        return []

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
        code = [L.VariableDecl("const double", access, expr)]
        return code

    def _expect_table(self, e, mt, tabledata, num_points, access):
        """These quantities refer to constant tables defined in ufc_geometry.h."""
        # TODO: Inject const static table here instead?
        return []

    reference_cell_volume = _expect_table
    reference_facet_volume = _expect_table
    reference_normal = _expect_table
    cell_facet_jacobian = _expect_table
    reference_cell_edge_vectors = _expect_table
    reference_facet_edge_vectors = _expect_table
    facet_orientation = _expect_table

    def _expect_physical_coords(self, e, mt, tabledata, num_points, access):
        """These quantities refer to coordinate_dofs"""
        # TODO: Generate more efficient inline code for Max/MinCell/FacetEdgeLength
        #       and CellDiameter here rather than lowering these quantities?
        return []

    cell_vertices = _expect_physical_coords
    cell_edge_vectors = _expect_physical_coords
    facet_edge_vectors = _expect_physical_coords

    def _expect_symbolic_lowering(self, e, mt, tabledata, num_points, access):
        """These quantities are expected to be replaced in symbolic preprocessing."""
        logging.exception("Expecting {0} to be replaced in symbolic preprocessing.".format(type(e)))

    facet_normal = _expect_symbolic_lowering
    cell_normal = _expect_symbolic_lowering
    jacobian_inverse = _expect_symbolic_lowering
    jacobian_determinant = _expect_symbolic_lowering
    facet_jacobian = _expect_symbolic_lowering
    facet_jacobian_inverse = _expect_symbolic_lowering
    facet_jacobian_determinant = _expect_symbolic_lowering
