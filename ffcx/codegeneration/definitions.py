# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCx/UFC specific variable definitions."""

import logging

from typing import List
import ufl
import ffcx.codegeneration.lnodes as L

logger = logging.getLogger("ffcx")


def create_nested_for_loops(indices: List[L.MultiIndex], body):
    ranges = [r for idx in indices for r in idx.sizes]
    indices = [idx.local_index(i) for idx in indices for i in range(len(idx.sizes))]
    depth = len(ranges)
    for i in reversed(range(depth)):
        body = L.ForRange(indices[i], 0, ranges[i], body=[body])
    return body


def create_quadrature_index(quadrature_rule, quadrature_index_symbol):
    """Create a multi index for the quadrature loop."""
    ranges = [0]
    name = quadrature_index_symbol.name
    if quadrature_rule:
        indices = [L.Symbol(name, dtype=L.DataType.INT)]
        ranges = [quadrature_rule.weights.size]
        if quadrature_rule.has_tensor_factors:
            dim = len(quadrature_rule.tensor_factors)
            ranges = [factor[1].size for factor in quadrature_rule.tensor_factors]
            indices = [L.Symbol(name + f"{i}", dtype=L.DataType.INT) for i in range(dim)]

    return L.MultiIndex(indices, ranges)


def create_dof_index(tabledata, dof_index_symbol):
    """Create a multi index for the coefficient dofs."""
    name = dof_index_symbol.name
    if tabledata.has_tensor_factorisation:
        dim = len(tabledata.tensor_factors)
        ranges = [factor.values.shape[-1] for factor in tabledata.tensor_factors]
        indices = [L.Symbol(f"{name}{i}", dtype=L.DataType.INT) for i in range(dim)]
    else:
        ranges = [tabledata.values.shape[-1]]
        indices = [L.Symbol(name, dtype=L.DataType.INT)]

    return L.MultiIndex(indices, ranges)


class FFCXBackendDefinitions(object):
    """FFCx specific code definitions."""

    def __init__(self, ir, symbols, options):
        # Store ir and options
        self.integral_type = ir.integral_type
        self.entitytype = ir.entitytype
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
        # For applying tensor product to coefficients, we need to know if the coefficient
        # has a tensor factorisation and if the quadrature rule has a tensor factorisation.
        # If both are true, we can apply the tensor product to the coefficient.

        apply_tensor_product = tabledata.has_tensor_factorisation and quadrature_rule.has_tensor_factors

        iq_symbol = self.symbols.quadrature_loop_index
        ic_symbol = self.symbols.coefficient_dof_sum_index

        iq = create_quadrature_index(quadrature_rule, iq_symbol)
        ic = create_dof_index(tabledata, ic_symbol)

        # Get properties of tables
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
        # FE = self.symbols.element_table(tabledata, self.entitytype, mt.restriction)
        FE = self.symbols.table_access(tabledata, self.entitytype, mt.restriction, iq, ic)

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
                pre_code += [L.ArrayDecl(dof_access.array, sizes=num_dofs)]
                pre_body = [L.Assign(dof_access, dof_access_map)]
                pre_code += [L.ForRange(ic, 0, num_dofs, pre_body)]
        else:
            # For bs == 1, the coefficient access has a stride of 1. e.g.: XXXXXX
            dof_access = self.symbols.coefficient_dof_access(mt.terminal, (ic.global_index) * bs + begin)

        body = [L.AssignAdd(access, dof_access * FE)]
        code += [L.VariableDecl(access, 0.0)]
        if apply_tensor_product:
            code += [create_nested_for_loops([ic], body)]
        else:
            code += [L.ForRange(ic, 0, num_dofs, body)]

        return pre_code, code

    def constant(self, t, mt, tabledata, quadrature_rule, access):
        # Constants are not defined within the kernel.
        # No definition is needed because access to them is directly
        # via symbol c[], i.e. as passed into the kernel.
        return [], []

    def _define_coordinate_dofs_lincomb(self, e, mt, tabledata, quadrature_rule, access):
        """Define x or J as a linear combination of coordinate dofs with given table data."""
        # Get properties of domain
        domain = ufl.domain.extract_unique_domain(mt.terminal)
        coordinate_element = domain.ufl_coordinate_element()
        num_scalar_dofs = coordinate_element.sub_element.dim

        num_dofs = tabledata.values.shape[3]
        begin = tabledata.offset

        assert num_scalar_dofs == num_dofs

        # Find table name
        ttype = tabledata.ttype

        assert ttype != "zeros"
        assert ttype != "ones"

        # Get access to element table
        ic_symbol = self.symbols.coefficient_dof_sum_index
        iq_symbol = self.symbols.quadrature_loop_index
        ic = create_dof_index(tabledata, ic_symbol)
        iq = create_quadrature_index(quadrature_rule, iq_symbol)
        FE = self.symbols.table_access(tabledata, self.entitytype, mt.restriction, iq, ic)

        dof_access = L.Symbol("coordinate_dofs", dtype=L.DataType.REAL)

        # coordinate dofs is always 3d
        dim = 3
        offset = 0
        if mt.restriction == "-":
            offset = num_scalar_dofs * dim

        code = []
        body = [L.AssignAdd(access, dof_access[ic.global_index * dim + begin + offset] * FE)]
        code += [L.VariableDecl(access, 0.0)]
        code += [create_nested_for_loops([ic], body)]

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
