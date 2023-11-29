# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCx/UFC specific symbol naming."""

import logging
import ufl
import ffcx.codegeneration.lnodes as L

logger = logging.getLogger("ffcx")


# TODO: Get restriction postfix from somewhere central
def ufcx_restriction_postfix(restriction):
    if restriction == "+":
        res = "_0"
    elif restriction == "-":
        res = "_1"
    else:
        res = ""
    return res


def format_mt_name(basename, mt):
    """Format variable name for modified terminal."""
    access = str(basename)

    # Add averaged state to name
    if mt.averaged is not None:
        avg = f"_a{mt.averaged}"
        access += avg

    # Format restriction
    res = ufcx_restriction_postfix(mt.restriction).replace("_", "_r")
    access += res

    # Format global derivatives
    if mt.global_derivatives:
        assert basename == "J"
        der = f"_deriv_{''.join(map(str, mt.global_derivatives))}"
        access += der

    # Format local derivatives
    if mt.local_derivatives:
        # Convert "listing" derivative multindex into "counting" representation
        gdim = ufl.domain.extract_unique_domain(mt.terminal).geometric_dimension()
        ld_counting = tuple(mt.local_derivatives.count(i) for i in range(gdim))
        der = f"_d{''.join(map(str, ld_counting))}"
        access += der

    # Add flattened component to name
    if mt.component:
        comp = f"_c{mt.flat_component}"
        access += comp

    return access


class FFCXBackendSymbols(object):
    """FFCx specific symbol definitions. Provides non-ufl symbols."""

    def __init__(self, coefficient_numbering, coefficient_offsets,
                 original_constant_offsets):
        self.coefficient_numbering = coefficient_numbering
        self.coefficient_offsets = coefficient_offsets

        self.original_constant_offsets = original_constant_offsets

        # Keep tabs on tables, so the symbols can be reused
        self.quadrature_weight_tables = {}
        self.element_tables = {}

        # Reusing a single symbol for all quadrature loops, assumed not to be nested.
        self.quadrature_loop_index = L.Symbol("iq", dtype=L.DataType.INT)

        # Symbols for the tabulate_tensor function arguments
        self.element_tensor = L.Symbol("A", dtype=L.DataType.SCALAR)
        self.coefficients = L.Symbol("w", dtype=L.DataType.SCALAR)
        self.constants = L.Symbol("c", dtype=L.DataType.SCALAR)
        self.coordinate_dofs = L.Symbol("coordinate_dofs", dtype=L.DataType.REAL)
        self.entity_local_index = L.Symbol("entity_local_index", dtype=L.DataType.INT)
        self.quadrature_permutation = L.Symbol("quadrature_permutation", dtype=L.DataType.INT)

        # Index for loops over coefficient dofs, assumed to never be used in two nested loops.
        self.coefficient_dof_sum_index = L.Symbol("ic", dtype=L.DataType.INT)

        # Table for chunk of custom quadrature weights (including cell measure scaling).
        self.custom_weights_table = L.Symbol("weights_chunk", dtype=L.DataType.REAL)

        # Table for chunk of custom quadrature points (physical coordinates).
        self.custom_points_table = L.Symbol("points_chunk", dtype=L.DataType.REAL)

    def entity(self, entitytype, restriction):
        """Entity index for lookup in element tables."""
        if entitytype == "cell":
            # Always 0 for cells (even with restriction)
            return L.LiteralInt(0)

        if entitytype == "facet":
            if restriction == "-":
                return self.entity_local_index[1]
            else:
                return self.entity_local_index[0]
        elif entitytype == "vertex":
            return self.entity_local_index[0]
        else:
            logging.exception(f"Unknown entitytype {entitytype}")

    def argument_loop_index(self, iarg):
        """Loop index for argument #iarg."""
        indices = ["i", "j", "k", "l"]
        return L.Symbol(indices[iarg], dtype=L.DataType.INT)

    def weights_table(self, quadrature_rule):
        """Table of quadrature weights."""
        key = f"weights_{quadrature_rule.id()}"
        if key not in self.quadrature_weight_tables:
            self.quadrature_weight_tables[key] = L.Symbol(f"weights_{quadrature_rule.id()}",
                                                          dtype=L.DataType.REAL)
        return self.quadrature_weight_tables[key]

    def points_table(self, quadrature_rule):
        """Table of quadrature points (points on the reference integration entity)."""
        return L.Symbol(f"points_{quadrature_rule.id()}", dtype=L.DataType.REAL)

    def x_component(self, mt):
        """Physical coordinate component."""
        return L.Symbol(format_mt_name("x", mt), dtype=L.DataType.REAL)

    def J_component(self, mt):
        """Jacobian component."""
        # FIXME: Add domain number!
        return L.Symbol(format_mt_name("J", mt), dtype=L.DataType.REAL)

    def domain_dof_access(self, dof, component, gdim, num_scalar_dofs, restriction):
        # FIXME: Add domain number or offset!
        offset = 0
        if restriction == "-":
            offset = num_scalar_dofs * 3
        return self.coordinate_dofs[3 * dof + component + offset]

    def domain_dofs_access(self, gdim, num_scalar_dofs, restriction):
        # FIXME: Add domain number or offset!
        return [
            self.domain_dof_access(dof, component, gdim, num_scalar_dofs, restriction)
            for dof in range(num_scalar_dofs) for component in range(gdim)
        ]

    def coefficient_dof_access(self, coefficient, dof_index):
        offset = self.coefficient_offsets[coefficient]
        w = self.coefficients
        return w[offset + dof_index]

    def coefficient_dof_access_blocked(self, coefficient: ufl.Coefficient, index,
                                       block_size, dof_offset):
        coeff_offset = self.coefficient_offsets[coefficient]
        w = self.coefficients
        _w = L.Symbol(f"_w_{coeff_offset}_{dof_offset}", dtype=L.DataType.SCALAR)
        unit_stride_access = _w[index]
        original_access = w[coeff_offset + index * block_size + dof_offset]
        return unit_stride_access, original_access

    def coefficient_value(self, mt):
        """Symbol for variable holding value or derivative component of coefficient."""
        c = self.coefficient_numbering[mt.terminal]
        return L.Symbol(format_mt_name("w%d" % (c, ), mt), dtype=L.DataType.SCALAR)

    def constant_index_access(self, constant, index):
        offset = self.original_constant_offsets[constant]
        c = self.constants
        return c[offset + index]

    def element_table(self, tabledata, entitytype, restriction):
        entity = self.entity(entitytype, restriction)

        if tabledata.is_uniform:
            entity = 0
        else:
            entity = self.entity(entitytype, restriction)

        if tabledata.is_piecewise:
            iq = 0
        else:
            iq = self.quadrature_loop_index

        if tabledata.is_permuted:
            qp = self.quadrature_permutation[0]
            if restriction == "-":
                qp = self.quadrature_permutation[1]
        else:
            qp = 0

        # Return direct access to element table, reusing symbol if possible
        if tabledata.name not in self.element_tables:
            self.element_tables[tabledata.name] = L.Symbol(tabledata.name,
                                                           dtype=L.DataType.REAL)
        return self.element_tables[tabledata.name][qp][entity][iq]
