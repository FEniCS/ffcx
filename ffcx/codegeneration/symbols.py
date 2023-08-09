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


class FFCXSymbols(object):
    """FFCx specific symbol definitions. Provides non-ufl symbols."""

    def __init__(self, ir):
        self.coefficient_numbering = ir.coefficient_numbering
        self.coefficient_offsets = ir.coefficient_offsets
        self.original_constant_offsets = ir.original_constant_offsets

    def element_tensor(self):
        """Symbol for the element tensor itself."""
        return L.Symbol("A")

    def entity(self, entitytype, restriction):
        """Entity index for lookup in element tables."""
        if entitytype == "cell":
            # Always 0 for cells (even with restriction)
            return L.LiteralInt(0)
        elif entitytype == "facet":
            postfix = "[0]"
            if restriction == "-":
                postfix = "[1]"
            return L.Symbol("entity_local_index" + postfix)
        elif entitytype == "vertex":
            return L.Symbol("entity_local_index[0]")
        else:
            logging.exception(f"Unknown entitytype {entitytype}")

    def argument_loop_index(self, iarg):
        """Loop index for argument #iarg."""
        indices = ["i", "j", "k", "l"]
        return L.Symbol(indices[iarg])

    def coefficient_dof_sum_index(self):
        """Index for loops over coefficient dofs, assumed to never be used in two nested loops."""
        return L.Symbol("ic")

    def quadrature_loop_index(self):
        """Reusing a single index name for all quadrature loops, assumed not to be nested."""
        return L.Symbol("iq")

    def quadrature_permutation(self, index):
        """Quadrature permutation, as input to the function."""
        return L.Symbol("quadrature_permutation")[index]

    def custom_weights_table(self):
        """Table for chunk of custom quadrature weights (including cell measure scaling)."""
        return L.Symbol("weights_chunk")

    def custom_points_table(self):
        """Table for chunk of custom quadrature points (physical coordinates)."""
        return L.Symbol("points_chunk")

    def weights_table(self, quadrature_rule):
        """Table of quadrature weights."""
        return L.Symbol(f"weights_{quadrature_rule.id()}")

    def points_table(self, quadrature_rule):
        """Table of quadrature points (points on the reference integration entity)."""
        return L.Symbol(f"points_{quadrature_rule.id()}")

    def x_component(self, mt):
        """Physical coordinate component."""
        return L.Symbol(format_mt_name("x", mt))

    def J_component(self, mt):
        """Jacobian component."""
        # FIXME: Add domain number!
        return L.Symbol(format_mt_name("J", mt))

    def domain_dof_access(self, dof, component, gdim, num_scalar_dofs, restriction):
        # FIXME: Add domain number or offset!
        offset = 0
        if restriction == "-":
            offset = num_scalar_dofs * 3
        vc = L.Symbol("coordinate_dofs")
        return vc[3 * dof + component + offset]

    def domain_dofs_access(self, gdim, num_scalar_dofs, restriction):
        # FIXME: Add domain number or offset!
        return [
            self.domain_dof_access(dof, component, gdim, num_scalar_dofs, restriction)
            for dof in range(num_scalar_dofs)
            for component in range(gdim)
        ]

    def coefficient_dof_access(self, coefficient, dof_index):
        offset = self.coefficient_offsets[coefficient]
        w = L.Symbol("w")
        return w[offset + dof_index]

    def coefficient_dof_access_blocked(
        self, coefficient: ufl.Coefficient, index, block_size, dof_offset
    ):
        coeff_offset = self.coefficient_offsets[coefficient]
        w = L.Symbol("w")
        _w = L.Symbol(f"_w_{coeff_offset}_{dof_offset}")
        unit_stride_access = _w[index]
        original_access = w[coeff_offset + index * block_size + dof_offset]
        return unit_stride_access, original_access

    def coefficient_value(self, mt):
        """Symbol for variable holding value or derivative component of coefficient."""
        c = self.coefficient_numbering[mt.terminal]
        return L.Symbol(format_mt_name("w%d" % (c,), mt))

    def constant_index_access(self, constant, index):
        offset = self.original_constant_offsets[constant]
        c = L.Symbol("c")
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
            iq = self.quadrature_loop_index()

        if tabledata.is_permuted:
            qp = self.quadrature_permutation(0)
            if restriction == "-":
                qp = self.quadrature_permutation(1)
        else:
            qp = 0

        # Return direct access to element table
        return L.Symbol(tabledata.name)[qp][entity][iq]
