# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCX/UFC specific symbol naming."""

import logging

logger = logging.getLogger(__name__)


# TODO: Get restriction postfix from somewhere central
def ufc_restriction_postfix(restriction):
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
    if mt.averaged:
        avg = "_a{0}".format(mt.averaged)
        access += avg

    # Format restriction
    res = ufc_restriction_postfix(mt.restriction).replace("_", "_r")
    access += res

    # Format local derivatives
    assert not mt.global_derivatives
    if mt.local_derivatives:
        der = "_d{0}".format(''.join(map(str, mt.local_derivatives)))
        access += der

    # Add flattened component to name
    if mt.component:
        comp = "_c{0}".format(mt.flat_component)
        access += comp

    return access


class FFCXBackendSymbols(object):
    """FFCX specific symbol definitions. Provides non-ufl symbols."""

    def __init__(self, language, coefficient_numbering, coefficient_offsets,
                 original_constant_offsets):
        self.L = language
        self.S = self.L.Symbol
        self.coefficient_numbering = coefficient_numbering
        self.coefficient_offsets = coefficient_offsets

        self.original_constant_offsets = original_constant_offsets

        # Used for padding variable names based on restriction
#        self.restriction_postfix = {r: ufc_restriction_postfix(r) for r in ("+", "-", None)}

        # TODO: Make this configurable for easy experimentation with dolfinx!
        # Coordinate dofs for each component are interleaved? Must match dolfinx.
        # True = XYZXYZXYZXYZ, False = XXXXYYYYZZZZ
        self.interleaved_components = True

    def element_tensor(self):
        """Symbol for the element tensor itself."""
        return self.S("A")

    def entity(self, entitytype, restriction):
        """Entity index for lookup in element tables."""
        if entitytype == "cell":
            # Always 0 for cells (even with restriction)
            return self.L.LiteralInt(0)
        elif entitytype == "facet":
            postfix = "[0]"
            if restriction == "-":
                postfix = "[1]"
            return self.S("facet" + postfix)
        elif entitytype == "vertex":
            return self.S("vertex[0]")
        else:
            logging.exception("Unknown entitytype {}".format(entitytype))

    def argument_loop_index(self, iarg):
        """Loop index for argument #iarg."""
        indices = ["i", "j", "k", "l"]
        return self.S(indices[iarg])

    def entity_reflection(self, L, i, cell_shape):
        """Returns the bool that says whether or not an entity has been reflected."""
        cell_info = self.S("cell_permutation")
        if cell_shape in ["triangle", "quadrilateral"]:
            num_faces = 0
            face_bitsize = 1
            assert i[0] == 1
        if cell_shape == "tetrahedron":
            num_faces = 4
            face_bitsize = 3
        if cell_shape == "hexahedron":
            num_faces = 4
            face_bitsize = 3
        if i[0] == 1:
            return L.NE(L.BitwiseAnd(cell_info, L.BitShiftL(1, face_bitsize * num_faces + i[1])), 0)
        elif i[0] == 2:
            return L.NE(L.BitwiseAnd(cell_info, L.BitShiftL(1, face_bitsize * i[1])), 0)
        return L.LiteralBool(False)

    def entity_rotations(self, L, i, cell_shape):
        """Returns the bool that says whether or not an entity has been reflected."""
        cell_info = self.S("cell_permutation")
        if cell_shape == "tetrahedron":
            face_bitsize = 3
        if cell_shape == "hexahedron":
            face_bitsize = 3
        if i[0] == 2:
            return L.BitwiseAnd(cell_info, L.BitShiftL(3, face_bitsize * i[1] + 1))
        return L.LiteralBool(False)

    def coefficient_dof_sum_index(self):
        """Index for loops over coefficient dofs, assumed to never be used in two nested loops."""
        return self.S("ic")

    def quadrature_loop_index(self):
        """Reusing a single index name for all quadrature loops, assumed not to be nested."""
        return self.S("iq")

    def quadrature_permutation(self, index):
        """Quadrature permutation, as input to the function."""
        return self.S("quadrature_permutation")[index]

    def num_custom_quadrature_points(self):
        """Number of quadrature points, argument to custom integrals."""
        return self.S("num_quadrature_points")

    def custom_quadrature_weights(self):
        """Quadrature weights including cell measure scaling, argument to custom integrals."""
        return self.S("quadrature_weights")

    def custom_quadrature_points(self):
        """Physical quadrature points, argument to custom integrals."""
        return self.S("quadrature_points")

    def custom_weights_table(self):
        """Table for chunk of custom quadrature weights (including cell measure scaling)."""
        return self.S("weights_chunk")

    def custom_points_table(self):
        """Table for chunk of custom quadrature points (physical coordinates)."""
        return self.S("points_chunk")

    def weights_table(self, num_points):
        """Table of quadrature weights."""
        return self.S("weights%d" % (num_points, ))

    def points_table(self, num_points):
        """Table of quadrature points (points on the reference integration entity)."""
        return self.S("points%d" % (num_points, ))

    def x_component(self, mt):
        """Physical coordinate component."""
        return self.S(format_mt_name("x", mt))

    def J_component(self, mt):
        """Jacobian component."""
        # FIXME: Add domain number!
        return self.S(format_mt_name("J", mt))

    def domain_dof_access(self, dof, component, gdim, num_scalar_dofs, restriction):
        # FIXME: Add domain number or offset!
        offset = 0
        if restriction == "-":
            offset = num_scalar_dofs * gdim
        vc = self.S("coordinate_dofs")
        if self.interleaved_components:
            return vc[gdim * dof + component + offset]
        else:
            return vc[num_scalar_dofs * component + dof + offset]

    def domain_dofs_access(self, gdim, num_scalar_dofs, restriction):
        # FIXME: Add domain number or offset!
        return [
            self.domain_dof_access(dof, component, gdim, num_scalar_dofs, restriction)
            for component in range(gdim) for dof in range(num_scalar_dofs)
        ]

    def coefficient_dof_access(self, coefficient, dof_number):
        # TODO: Add domain number?
        offset = self.coefficient_offsets[coefficient]
        w = self.S("w")
        return w[offset + dof_number]

    def coefficient_value(self, mt):
        """Symbol for variable holding value or derivative component of coefficient."""

        c = self.coefficient_numbering[mt.terminal]
        return self.S(format_mt_name("w%d" % (c, ), mt))

    def constant_index_access(self, constant, index):
        offset = self.original_constant_offsets[constant]
        c = self.S("c")

        return c[offset + index]

    def named_table(self, name):
        return self.S(name)

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
        return self.named_table(tabledata.name)[qp][entity][iq]

    def expr_component_index(self):
        """Symbol for indexing the expression's ufl shape."""
        return self.S("ec")
