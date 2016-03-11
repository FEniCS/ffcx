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


from ufl.cell import num_cell_entities

def num_coordinate_component_dofs(coordinate_element):
    """Get the number of dofs for a coordinate component for this degree.

    This is a local hack that works for Lagrange 1-3, better
    would be to get this passed by ffc from fiat through the ir.
    The table data is to messy to figure out a clean design for that quickly.
    """
    degree = coordinate_element.degree()
    cell = coordinate_element.cell()
    tdim = cell.topological_dimension()
    cellname = cell.cellname()
    d = 0
    for entity_dim in range(tdim+1):
        # n = dofs per cell entity
        if entity_dim == 0:
            n = 1
        elif entity_dim == 1:
            n = degree - 1
        elif entity_dim == 2:
            n = (degree - 2)*(degree - 1) // 2
        elif entity_dim == 3:
            n = (degree - 3)*(degree - 2)*(degree - 1) // 6
        else:
            error("Entity dimension out of range")
        # Accumulate
        num_entities = num_cell_entities[cellname][entity_dim]
        d += num_entities * n
    return d


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
        return []

    def expr(self, t, mt, tabledata, access):
        error("Unhandled type {0}".format(type(t)))

    # === Generate code definitions ===

    def quadrature_weight(self, e, mt, tabledata, access):
        return []

    def constant_value(self, e, mt, tabledata, access):
        return []

    def argument(self, t, mt, tabledata, access):
        return []

    def coefficient(self, t, mt, tabledata, access):
        L = self.language

        # For a constant coefficient we reference the dofs directly, so no definition needed
        if is_cellwise_constant(mt.terminal):
            return []

        # No need to store basis function value in its own variable,
        # just get table value directly
        uname, begin, end = tabledata
        uname = L.Symbol(uname)

        # Empty loop needs to be skipped as zero tables may not be generated
        # FIXME: assert begin < end instead, and remove at earlier
        #        stage so dependent code can also be removed
        if begin >= end:
            return []

        # Get various symbols
        entity = self.symbols.entity(self.ir["entitytype"], mt.restriction)
        iq = self.symbols.quadrature_loop_index()
        idof = self.symbols.coefficient_dof_sum_index()
        dof_access = self.symbols.coefficient_dof_access(mt.terminal, idof)
        table_access = uname[entity][iq][idof - begin]

        # Loop to accumulate linear combination of dofs and tables
        code = [
            L.VariableDecl("double", access, 0.0),
            L.ForRange(idof, begin, end,
                       body=[L.AssignAdd(access, dof_access * table_access)])
            ]
        return code

    def _define_coordinate_dofs_lincomb(self, e, mt, tabledata, access):
        "Define something (x or J) linear combination of coordinate dofs with given table data."

        L = self.language

        # Get properties of domain
        domain = mt.terminal.ufl_domain()
        tdim = domain.topological_dimension()
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()
        degree = coordinate_element.degree()
        num_scalar_dofs = num_coordinate_component_dofs(coordinate_element)

        # Reference coordinates are known, no coordinate field, so we compute
        # this component as linear combination of coordinate_dofs "dofs" and table

        uname, begin, end = tabledata
        uname = L.Symbol(uname)
        assert end - begin <= num_scalar_dofs

        entity = self.symbols.entity(self.ir["entitytype"], mt.restriction)

        if is_cellwise_constant(mt.expr):
            iq = 0
        else:
            iq = self.symbols.quadrature_loop_index()

        if 0:  # FIXME: Make an option to test
            # Generated loop version:
            coefficient_dof = self.symbols.coefficient_dof_sum_index()
            dof_access = self.symbols.domain_dof_access(coefficient_dof, mt.flat_component,
                                                        gdim, num_scalar_dofs,
                                                        mt.restriction, self.interleaved_components)
            table_access = uname[entity][iq][coefficient_dof]

            # Loop to accumulate linear combination of dofs and tables
            code = [
                L.VariableDecl("double", access, 0.0),
                L.ForRange(coefficient_dof, begin, end,
                           body=[L.AssignAdd(access, dof_access * table_access)])
                ]
        else:
            # Inlined version (we know this is bounded by a small number)
            dof_access = self.symbols.domain_dofs_access(gdim, num_scalar_dofs,
                                                         mt.restriction,
                                                         self.interleaved_components)

            value = L.Sum([dof_access[idof] * uname[entity][iq][idof - begin]
                           for idof in range(begin, end)])

            # Inlined loop to accumulate linear combination of dofs and tables
            code = [
                L.VariableDecl("const double", access, value)
                ]

        return code

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
            return self._define_coordinate_dofs_lincomb(e, mt, tabledata, access)

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
        """Return definition code for the Jacobian of x(X).

        J = sum_k xdof_k grad_X xphi_k(X)
        """
        if self.physical_coordinates_known:
            return []
        else:
            return self._define_coordinate_dofs_lincomb(e, mt, tabledata, access)

    def cell_orientation(self, e, mt, tabledata, access):
        # Would be nicer if cell_orientation was a double variable input,
        # but this is how dolfin/ufc/ffc currently passes this information.
        # 0 means up and gives +1.0, 1 means down and gives -1.0.
        L = self.language
        co = "cell_orientation" + ufc_restriction_postfix(mt.restriction)
        expr = L.VerbatimExpr("(" + co + " == 1) ? -1.0: +1.0;")
        code = [
            L.VariableDecl("const double", access, expr)
            ]
        return code

    def _expect_table(self, e, mt, tabledata, access):
        "These quantities refer to constant tables defined in ufc_geometry.h."
        # TODO: Inject const static table here instead?
        return []
    reference_cell_volume = _expect_table
    reference_facet_volume = _expect_table
    reference_normal = _expect_table
    cell_facet_jacobian = _expect_table
    cell_edge_vectors = _expect_table
    facet_edge_vectors = _expect_table
    facet_orientation = _expect_table

    def _expect_symbolic_lowering(self, e, mt, tabledata, access):
        "These quantities are expected to be replaced in symbolic preprocessing."
        error("Expecting {0} to be replaced in symbolic preprocessing.".format(type(e)))
    facet_normal = _expect_symbolic_lowering
    cell_normal = _expect_symbolic_lowering
    jacobian_inverse = _expect_symbolic_lowering
    jacobian_determinant = _expect_symbolic_lowering
    facet_jacobian = _expect_symbolic_lowering
    facet_jacobian_inverse = _expect_symbolic_lowering
    facet_jacobian_determinant = _expect_symbolic_lowering
