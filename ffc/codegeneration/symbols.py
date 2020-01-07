# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFC/UFC specific symbol naming."""

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


class FFCBackendSymbols(object):
    """FFC specific symbol definitions. Provides non-ufl symbols."""

    def __init__(self, language, coefficient_numbering, coefficient_offsets,
                 original_constant_offsets):
        self.L = language
        self.S = self.L.Symbol
        self.coefficient_numbering = coefficient_numbering
        self.coefficient_offsets = coefficient_offsets

        self.original_constant_offsets = original_constant_offsets

        # Used for padding variable names based on restriction
#        self.restriction_postfix = {r: ufc_restriction_postfix(r) for r in ("+", "-", None)}

        # TODO: Make this configurable for easy experimentation with dolfin!
        # Coordinate dofs for each component are interleaved? Must match dolfin.
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

    def cell_orientation_argument(self, restriction):
        """Cell orientation argument in ufc. Not same as cell orientation in generated code."""
        postfix = "[0]"
        if restriction == "-":
            postfix = "[1]"
        return self.S("cell_orientation" + postfix)

    def cell_orientation_internal(self, restriction):
        """Internal value for cell orientation in generated code."""
        return self.S("co" + ufc_restriction_postfix(restriction))

    def argument_loop_index(self, iarg):
        """Loop index for argument #iarg."""
        indices = ["i", "j", "k", "l"]
        return self.S(indices[iarg])

    def coefficient_dof_sum_index(self):
        """Index for loops over coefficient dofs, assumed to never be used in two nested loops."""
        return self.S("ic")

    def quadrature_loop_index(self):
        """Reusing a single index name for all quadrature loops, assumed not to be nested."""
        return self.S("iq")

    def permuted_quadrature_loop_index(self, offset, num_points, celltype):
        """Reusing a single index name for all quadrature loops, assumed not to be nested."""
        iq = self.quadrature_loop_index()
        perm = self.quadrature_permutation()
        if celltype == "line":
            return self.L.Conditional(self.L.EQ(perm[offset], 1), num_points - 1 - iq, iq)
        if celltype == "square" or celltype == "triangle":
            return self.rotation_of_point(iq, num_points, (perm[offset]-2)//2, (perm[offset]-2)%2, celltype)
        return iq

    def rotation_of_point(self, iq, num_points, num_rotations, num_reflections, celltype):
        """Maybe the worst code I've ever written. Please replace me"""
        if celltype == "square":
            rotation, reflection = self.square_rotation_of_point(iq, num_points, num_rotations, num_reflections)
        if celltype == "triangle":
            rotation, reflection = self.triangle_rotation_of_point(iq, num_points, num_rotations, num_reflections)

        assert len(rotation) == len(reflection) == num_points

        permed = ""
        ls = list(range(num_points))
        for i in range(num_points):
            iq_p = ""
            for j in range(num_points):
                if_ref = self.L.Conditional(self.L.EQ(num_reflections, 1), reflection[ls[j]], ls[j])
                if iq_p == "":
                    iq_p = if_ref
                else:
                    iq_p = self.L.Conditional(self.L.EQ(iq, j), if_ref, iq_p)
            if permed == "":
                permed = iq_p
            else:
                permed = self.L.Conditional(self.L.EQ(num_rotations, i), iq_p, permed)
            ls = [rotation[a] for a in ls]
            #ls = [ls[a] for a in rotation]
        return permed

    def square_rotation_of_point(self, iq, num_points, num_rotations, num_reflections):
        """Maybe the worst code I've ever written. Please replace me"""
        #### TODO: look at generated code to debug this...
        side = 0
        while side ** 2 < num_points:
            side += 1
        assert num_points == side ** 2
        rotation = [side - 1 - j + side * i for j in range(side) for i in range(side)]
        reflection = [j + side * i for j in range(side) for i in range(side)]
        return rotation, reflection

    def triangle_rotation_of_point(self, iq, num_points, num_rotations, num_reflections):
        """Maybe the worst code I've ever written. Please replace me"""
        side = 0
        while side * (side + 1) //  2 < num_points:
            side += 1
        assert num_points == side * (side + 1) // 2
        rotation = [side - j - 1 + i * side - i * (i + 1) // 2 for j in range(side) for i in range(side-j)]
        reflection = [j + i * (side + 1) - i * (i + 1) // 2 for j in range(side) for i in range(side-j)]
        return rotation, reflection

    def quadrature_permutation(self):
        """Quadrature permutation, as input to the function."""
        return self.S("quadrature_permutation")

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

    def element_table(self, tabledata, entitytype, restriction):
        if tabledata.is_uniform:
            entity = 0
        else:
            entity = self.entity(entitytype, restriction)

        if tabledata.is_piecewise:
            iq = 0
        else:
            iq = self.quadrature_loop_index()

        # Return direct access to element table
        return self.S(tabledata.name)[entity][iq]

    def element_table_flip(self, tabledata, entitytype, restriction, num_points, celltype):
        if tabledata.is_uniform:
            entity = 0
        else:
            entity = self.entity(entitytype, restriction)

        if tabledata.is_piecewise:
            iq = 0
        else:
            if restriction == "+":
                iq = self.permuted_quadrature_loop_index(0, num_points, celltype)
            else:
                iq = self.permuted_quadrature_loop_index(1, num_points, celltype)

        # Return direct access to element table
        return self.S(tabledata.name)[entity][iq]

    def expr_component_index(self):
        """Symbol for indexing the expression's ufl shape."""
        return self.S("ec")
