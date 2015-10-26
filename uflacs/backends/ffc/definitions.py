# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""FFC specific definitions."""

from six.moves import xrange as range

from ufl.corealg.multifunction import MultiFunction
from ufl.checks import is_cellwise_constant

from ffc.log import error
from ffc.log import ffc_assert

from uflacs.backends.ffc.common import FFCBackendSymbols
# FIXME: Move these to FFCBackendSymbols
from uflacs.backends.ffc.common import format_entity_name, ufc_restriction_postfix


class FFCDefinitionsBackend(MultiFunction):
    """FFC specific code definitions."""

    def __init__(self, ir, language, parameters):
        MultiFunction.__init__(self)

        # Store ir and parameters
        self.ir = ir
        self.language = language
        self.parameters = parameters

        # FIXME: Make this configurable for easy experimentation with dolfin!
        # Coordinate dofs for each component are interleaved? Must match dolfin.
        self.interleaved_components = True # parameters["interleaved_coordinate_component_dofs"]

        # Configure definitions behaviour
        if self.ir["integral_type"] in ("custom", "vertex"):
            self.physical_coordinates_known = True
        else:
            self.physical_coordinates_known = False

        # Need this for custom integrals
        #classname = make_classname(prefix, "finite_element", ir["element_numbers"][ufl_element])

        coefficient_numbering = ir["uflacs"]["coefficient_numbering"]
        self.symbols = FFCBackendSymbols(self.language, coefficient_numbering)

    def get_includes(self):
        "Return include statements to insert at top of file."
        includes = []
        return includes

    def initial(self):
        "Return code inserted at beginning of kernel."
        code = []
        return code

    def expr(self, t, mt, tabledata, access):
        error("Unhandled type {0}".format(type(t)))

    # === Generate code definitions ===

    def constant_value(self, e, mt, tabledata, access):
        return None

    def argument(self, t, mt, tabledata, access):
        code = []
        return code

    def coefficient(self, t, mt, tabledata, access):
        L = self.language

        code = []
        if is_cellwise_constant(mt.terminal):
            # For a constant coefficient we reference the dofs directly, so no definition needed
            pass
        else:
            # No need to store basis function value in its own variable,
            # just get table value directly
            code += [L.VariableDecl("double", access, 0.0)]
            uname, begin, end = tabledata
            entity = format_entity_name(self.ir["entitytype"], mt.restriction)

            # Empty loop needs to be skipped as zero tables may not be generated
            # FIXME: assert begin < end instead, and remove at earlier
            #        stage so dependent code can also be removed
            if begin >= end:
                return code

            iq = self.symbols.quadrature_loop_index()
            idof = self.symbols.coefficient_dof_sum_index()

            dof = L.Sub(idof, begin)
            table_access = L.ArrayAccess(uname, (entity, iq, dof))

            dof_access = self.symbols.coefficient_dof_access(mt.terminal, idof)

            prod = L.Mul(dof_access, table_access)
            body = [L.AssignAdd(access, prod)]

            # Loop to accumulate linear combination of dofs and tables
            code += [L.ForRange(idof, begin, end, body=body)]

        return code

    def quadrature_weight(self, e, mt, tabledata, access):
        return []

    def spatial_coordinate(self, e, mt, tabledata, access):
        """Return definition code for the physical spatial coordinates.

        If physical coordinates are given:
          No definition needed.

        If reference coordinates are given:
          x = sum_k xdof_k xphi_k(X)

        If reference facet coordinates are given:
          x = sum_k xdof_k xphi_k(Xf)
        """
        if self.physical_coordinates_known:
            return []
        else:
            return self._define_spatial_coordinate(e, mt, tabledata, access)

    def _define_spatial_coordinate(self, e, mt, tabledata, access):
        """Return definition code for x(X).

        x = sum_k xdof_k xphi_k(X)
        """
        # FIXME: Share code with _define_jacobian
        # access here is e.g. x0, component 0 of x

        L = self.language
        code = []

        coordinate_element = mt.terminal.ufl_domain().ufl_coordinate_element()
        degree = coordinate_element.degree()

        # Reference coordinates are known, no coordinate field, so we compute
        # this component as linear combination of coordinate_dofs "dofs" and table

        cell = mt.terminal.ufl_domain().ufl_cell()
        gdim = cell.geometric_dimension()

        uname, begin, end = tabledata

        if degree == 1:
            num_vertices = cell.num_vertices()
            # FIXME: This assertion fails if we use num_scalar_dofs=end-begin
            #ffc_assert(num_scalar_dofs <= num_vertices,
            #           "For a linear element, should have dofs(%d) == vertices(%d). (Note: comp=%d, begin=%d, end=%d)." % (
            #               num_scalar_dofs, num_vertices, mt.flat_component, begin, end))
            num_scalar_dofs = num_vertices
        else:
            num_scalar_dofs = end - begin

        entity = format_entity_name(self.ir["entitytype"], mt.restriction)
        coefficient_dof = self.symbols.coefficient_dof_sum_index()
        if coordinate_element.degree() > 0:
            iq = self.symbols.quadrature_loop_index()
        else:
            iq = 0

        if 0:
            # Generated loop version:
            table_access = L.ArrayAccess(uname, (entity, iq, coefficient_dof))
            dof_access = self.symbols.domain_dof_access(coefficient_dof, mt.flat_component,
                                                        gdim, num_scalar_dofs,
                                                        mt.restriction, self.interleaved_components)
            prod = L.Mul(dof_access, table_access)
            accumulate = [L.AssignAdd(access, prod)]

            # Loop to accumulate linear combination of dofs and tables
            code += [L.VariableDecl("double", access, 0.0)]
            code += [L.ForRange(coefficient_dof, begin, end, body=accumulate)]
        else:
            # Inlined version (we know this is bounded by a small number)
            dof_access = self.symbols.domain_dofs_access(gdim, num_scalar_dofs, mt.restriction, self.interleaved_components)
            prods = []
            for idof in range(begin, end):
                ind = (entity, iq, L.Sub(idof, begin))
                table_access = L.ArrayAccess(uname, ind)
                prods += [L.Mul(dof_access[idof], table_access)]
                # TODO: Shorter notation possible here and elsewhere if symbols are L.Symbols and not strings:
                #prods += [dof_access[idof] * uname[entity, iq, idof - begin]]

            # Inlined loop to accumulate linear combination of dofs and tables
            code += [L.VariableDecl("const double", access, L.Sum(prods))]

        return code

    def cell_coordinate(self, e, mt, tabledata, access):
        """Return definition code for the reference spatial coordinates.

        If reference coordinates are given:
          No definition needed.

        If physical coordinates are given and domain is affine:
          X = K*(x-x0)
        This is inserted symbolically.

        If physical coordinates are given and domain is non- affine:
          Not currently supported.
        """
        return []

    def jacobian(self, e, mt, tabledata, access):
        if self.physical_coordinates_known:
            return []
        else:
            return self._define_jacobian(e, mt, tabledata, access)

    def _define_jacobian(self, e, mt, tabledata, access):
        """Return definition code for the Jacobian of x(X).

        J = sum_k xdof_k grad_X xphi_k(X)
        """
        # FIXME: Share code with _define_spatial_coordinate
        # access here is e.g. J_0, component 0 of J

        L = self.language
        code = []

        coordinate_element = mt.terminal.ufl_domain().ufl_coordinate_element()
        degree = coordinate_element.degree()

        # Reference coordinates are known, no coordinate field, so we compute
        # this component as linear combination of coordinate_dofs "dofs" and table

        cell = mt.terminal.ufl_domain().ufl_cell()
        gdim = cell.geometric_dimension()

        uname, begin, end = tabledata

        if degree == 1:
            num_vertices = cell.num_vertices()
            # FIXME: This assertion fails if we use num_scalar_dofs=end-begin
            #ffc_assert(num_scalar_dofs <= num_vertices,
            #           "For a linear element, should have dofs(%d) == vertices(%d). (Note: comp=%d, begin=%d, end=%d)." % (
            #               num_scalar_dofs, num_vertices, mt.flat_component, begin, end))
            num_scalar_dofs = num_vertices
        else:
            num_scalar_dofs = end - begin

        entity = format_entity_name(self.ir["entitytype"], mt.restriction)
        coefficient_dof = self.symbols.coefficient_dof_sum_index()
        if degree > 1:
            iq = self.symbols.quadrature_loop_index()
        else:
            iq = 0

        if 0:
            # Generated loop version:
            table_access = L.ArrayAccess(uname, (entity, iq, coefficient_dof))
            dof_access = self.symbols.domain_dof_access(coefficient_dof, mt.flat_component,
                                                        gdim, num_scalar_dofs,
                                                        mt.restriction, self.interleaved_components)
            prod = L.Mul(dof_access, table_access)
            accumulate = L.AssignAdd(access, prod)

            # Loop to accumulate linear combination of dofs and tables
            code += [L.VariableDecl("double", access, 0.0)]
            code += [L.ForRange(coefficient_dof, begin, end, body=accumulate)]
        else:
            # Inlined version:
            prods = []
            dof_access = self.symbols.domain_dofs_access(gdim, num_scalar_dofs, mt.restriction, self.interleaved_components)
            for idof in range(begin, end):
                ind = (entity, iq, L.Sub(idof, begin))
                table_access = L.ArrayAccess(uname, ind)
                prods += [L.Mul(dof_access[idof], table_access)]

            # Inlined loop to accumulate linear combination of dofs and tables
            code += [L.VariableDecl("const double", access, L.Sum(prods))]


        return code

    def reference_cell_volume(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        return []

    def reference_facet_volume(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        return []

    def reference_normal(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        return []

    def cell_facet_jacobian(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        return []

    def cell_edge_vectors(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        return []

    def facet_edge_vectors(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        return []

    def cell_orientation(self, e, mt, tabledata, access):
        # Would be nicer if cell_orientation was a double variable input,
        # but this is how dolfin/ufc/ffc currently passes this information.
        # 0 means up and gives +1.0, 1 means down and gives -1.0.
        L = self.language
        co = "cell_orientation" + ufc_restriction_postfix(mt.restriction)
        expr = L.VerbatimExpr("(" + co + " == 1) ? -1.0: +1.0;")
        return [L.VariableDecl("const double", access, expr)]

    def facet_orientation(self, e, mt, tabledata, access):
        # Constant table defined in ufc_geometry.h
        return []

    def _expect_symbolic_lowering(self, e, mt, tabledata, access):
        error("Expecting {0} to be replaced in symbolic preprocessing.".format(type(e)))
    facet_normal = _expect_symbolic_lowering
    cell_normal = _expect_symbolic_lowering
    jacobian_inverse = _expect_symbolic_lowering
    jacobian_determinant = _expect_symbolic_lowering
    facet_jacobian = _expect_symbolic_lowering
    facet_jacobian_inverse = _expect_symbolic_lowering
    facet_jacobian_determinant = _expect_symbolic_lowering
