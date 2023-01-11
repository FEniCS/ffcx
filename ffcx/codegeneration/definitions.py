# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCx/UFC specific variable definitions."""

import logging

import ufl
from ffcx.element_interface import create_element
from ffcx.naming import scalar_to_value_type

logger = logging.getLogger("ffcx")


class FFCXBackendDefinitions(object):
    """FFCx specific code definitions."""

    def __init__(self, ir, language, symbols, options):
        # Store ir and options
        self.integral_type = ir.integral_type
        self.entitytype = ir.entitytype
        self.language = language
        self.symbols = symbols
        self.options = options

        self.ir = ir

        # Lookup table for handler to call when the "get" method (below) is
        # called, depending on the first argument type.
        self.call_lookup = {ufl.coefficient.Coefficient: self.coefficient,
                            ufl.constant.Constant: self.constant,
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
                            ufl.geometry.CellOrientation: self._pass,
                            ufl.geometry.FacetOrientation: self._expect_table,
                            ufl.geometry.SpatialCoordinate: self.spatial_coordinate}

    def get(self, t, mt, tabledata, quadrature_rule, access):
        # Call appropriate handler, depending on the type of t
        ttype = type(t)
        handler = self.call_lookup.get(ttype, False)

        if not handler:
            # Look for parent class types instead
            for k in self.call_lookup.keys():
                if isinstance(t, k):
                    handler = self.call_lookup[k]
                    break

        if handler:
            return handler(t, mt, tabledata, quadrature_rule, access)
        else:
            raise RuntimeError("Not handled: %s", ttype)

    def coefficient(self, t, mt, tabledata, quadrature_rule, access):
        """Return definition code for coefficients."""
        L = self.language

        ttype = tabledata.ttype
        num_dofs = tabledata.values.shape[3]
        bs = tabledata.block_size
        begin = tabledata.offset
        end = begin + bs * (num_dofs - 1) + 1

        if ttype == "zeros":
            logging.debug("Not expecting zero coefficients to get this far.")
            return [], []

        # For a constant coefficient we reference the dofs directly, so no definition needed
        if ttype == "ones" and end - begin == 1:
            return [], []

        assert begin < end

        # Get access to element table
        FE = self.symbols.element_table(tabledata, self.entitytype, mt.restriction)
        ic = self.symbols.coefficient_dof_sum_index()

        code = []
        pre_code = []

        if bs > 1 and not tabledata.is_piecewise:
            # For bs > 1, the coefficient access has a stride of bs. e.g.: XYZXYZXYZ
            # When memory access patterns are non-sequential, the number of cache misses increases.
            # In turn, it results in noticeably reduced performance.
            # In this case, we create temp arrays outside the quadrature to store the coefficients and
            # have a sequential access pattern.
            dof_access, dof_access_map = self.symbols.coefficient_dof_access_blocked(mt.terminal, ic, bs, begin)

            # If a map is necessary from stride 1 to bs, the code must be added before the quadrature loop.
            if dof_access_map:
                pre_code += [L.ArrayDecl(self.options["scalar_type"], dof_access.array, num_dofs)]
                pre_body = L.Assign(dof_access, dof_access_map)
                pre_code += [L.ForRange(ic, 0, num_dofs, pre_body)]
        else:
            dof_access = self.symbols.coefficient_dof_access(mt.terminal, ic * bs + begin)

        body = [L.AssignAdd(access, dof_access * FE[ic])]
        code += [L.VariableDecl(self.options["scalar_type"], access, 0.0)]
        code += [L.ForRange(ic, 0, num_dofs, body)]

        return pre_code, code

    def constant(self, t, mt, tabledata, quadrature_rule, access):
        # Constants are not defined within the kernel.
        # No definition is needed because access to them is directly
        # via symbol c[], i.e. as passed into the kernel.
        return [], []

    def _define_coordinate_dofs_lincomb(self, e, mt, tabledata, quadrature_rule, access):
        """Define x or J as a linear combination of coordinate dofs with given table data."""
        L = self.language

        # Get properties of domain
        domain = ufl.domain.extract_unique_domain(mt.terminal)
        coordinate_element = domain.ufl_coordinate_element()
        num_scalar_dofs = create_element(coordinate_element).sub_element.dim

        num_dofs = tabledata.values.shape[3]
        begin = tabledata.offset

        assert num_scalar_dofs == num_dofs

        # Find table name
        ttype = tabledata.ttype

        assert ttype != "zeros"
        assert ttype != "ones"

        # Get access to element table
        FE = self.symbols.element_table(tabledata, self.entitytype, mt.restriction)
        ic = self.symbols.coefficient_dof_sum_index()
        dof_access = self.symbols.S("coordinate_dofs")

        # coordinate dofs is always 3d
        dim = 3
        offset = 0
        if mt.restriction == "-":
            offset = num_scalar_dofs * dim

        value_type = scalar_to_value_type(self.options["scalar_type"])

        code = []
        body = [L.AssignAdd(access, dof_access[ic * dim + begin + offset] * FE[ic])]
        code += [L.VariableDecl(f"{value_type}", access, 0.0)]
        code += [L.ForRange(ic, 0, num_scalar_dofs, body)]

        return [], code

    def spatial_coordinate(self, e, mt, tabledata, quadrature_rule, access):
        """Return definition code for the physical spatial coordinates.

        If physical coordinates are given:
          No definition needed.

        If reference coordinates are given:
          x = sum_k xdof_k xphi_k(X)

        If reference facet coordinates are given:
          x = sum_k xdof_k xphi_k(Xf)
        """
        if self.integral_type in ufl.custom_integral_types:
            # FIXME: Jacobian may need adjustment for custom_integral_types
            if mt.local_derivatives:
                logging.exception("FIXME: Jacobian in custom integrals is not implemented.")
            return []
        else:
            return self._define_coordinate_dofs_lincomb(e, mt, tabledata, quadrature_rule, access)

    def jacobian(self, e, mt, tabledata, quadrature_rule, access):
        """Return definition code for the Jacobian of x(X)."""
        return self._define_coordinate_dofs_lincomb(e, mt, tabledata, quadrature_rule, access)

    def _expect_table(self, e, mt, tabledata, quadrature_rule, access):
        """Return quantities referring to constant tables defined in the generated code."""
        # TODO: Inject const static table here instead?
        return [], []

    def _expect_physical_coords(self, e, mt, tabledata, quadrature_rule, access):
        """Return quantities referring to coordinate_dofs."""
        # TODO: Generate more efficient inline code for Max/MinCell/FacetEdgeLength
        #       and CellDiameter here rather than lowering these quantities?
        return [], []

    def _pass(self, *args, **kwargs):
        """Return nothing."""
        return [], []
