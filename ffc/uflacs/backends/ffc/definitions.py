# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
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

"""FFC/UFC specific variable definitions."""

from ufl.corealg.multifunction import MultiFunction
from ufl.measure import custom_integral_types

from ffc.log import error, warning

from ffc.uflacs.backends.ffc.symbols import FFCBackendSymbols
from ffc.uflacs.backends.ffc.common import num_coordinate_component_dofs


class FFCBackendDefinitions(MultiFunction):
    """FFC specific code definitions."""
    def __init__(self, ir, language, symbols, parameters):
        MultiFunction.__init__(self)

        # Store ir and parameters
        self.integral_type = ir["integral_type"]
        self.entitytype = ir["entitytype"]
        self.language = language
        self.symbols = symbols
        self.parameters = parameters


    # === Generate code to define variables for ufl types ===

    def expr(self, t, mt, tabledata, num_points, access):
        error("Unhandled type {0}".format(type(t)))


    #def quadrature_weight(self, e, mt, tabledata, num_points, access):
    #    "Quadrature weights are precomputed and need no code."
    #    return []


    def constant_value(self, e, mt, tabledata, num_points, access):
        "Constants simply use literals in the target language."
        return []


    def argument(self, t, mt, tabledata, num_points, access):
        "Arguments are accessed through element tables."
        return []


    def coefficient(self, t, mt, tabledata, num_points, access):
        "Return definition code for coefficients."
        L = self.language

        ttype = tabledata.ttype
        begin, end = tabledata.dofrange

        #fe_classname = ir["classnames"]["finite_element"][t.ufl_element()]

        if ttype == "zeros":
            debug("Not expecting zero coefficients to get this far.")
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
        #unroll = True
        if unroll:
            # TODO: Could also use a generated constant dofmap here like in block code
            # Unrolled loop to accumulate linear combination of dofs and tables
            values = [self.symbols.coefficient_dof_access(mt.terminal, idof) * FE[i]
                      for i, idof in enumerate(tabledata.dofmap)]
            value = L.Sum(values)
            code = [
                L.VariableDecl("const double", access, value)
                ]
        else:
            # Loop to accumulate linear combination of dofs and tables
            ic = self.symbols.coefficient_dof_sum_index()
            dof_access = self.symbols.coefficient_dof_access(mt.terminal, ic + begin)
            code = [
                L.VariableDecl("double", access, 0.0),
                L.ForRange(ic, 0, end - begin,
                           body=[L.AssignAdd(access, dof_access * FE[ic])])
                ]
        return code


    def _define_coordinate_dofs_lincomb(self, e, mt, tabledata, num_points, access):
        "Define x or J as a linear combination of coordinate dofs with given table data."
        L = self.language

        # Get properties of domain
        domain = mt.terminal.ufl_domain()
        #tdim = domain.topological_dimension()
        gdim = domain.geometric_dimension()
        coordinate_element = domain.ufl_coordinate_element()
        #degree = coordinate_element.degree()
        num_scalar_dofs = num_coordinate_component_dofs(coordinate_element)

        # Reference coordinates are known, no coordinate field, so we compute
        # this component as linear combination of coordinate_dofs "dofs" and table

        # Find table name and dof range it corresponds to
        ttype = tabledata.ttype
        begin, end = tabledata.dofrange

        assert end - begin <= num_scalar_dofs
        assert ttype != "zeros"
        assert ttype != "quadrature"
        #xfe_classname = ir["classnames"]["finite_element"][coordinate_element]
        #sfe_classname = ir["classnames"]["finite_element"][coordinate_element.sub_elements()[0]]

        # Get access to element table
        FE = self.symbols.element_table(tabledata, self.entitytype, mt.restriction)

        inline = True

        if ttype == "zeros":
            # Not sure if this will ever happen
            debug("Not expecting zeros for %s." % (e._ufl_class_.__name__,))
            code = [
                L.VariableDecl("const double", access, L.LiteralFloat(0.0))
                ]
        elif ttype == "ones":
            # Not sure if this will ever happen
            debug("Not expecting ones for %s." % (e._ufl_class_.__name__,))
            # Inlined version (we know this is bounded by a small number)
            dof_access = self.symbols.domain_dofs_access(gdim, num_scalar_dofs,
                                                         mt.restriction)
            values = [dof_access[idof] for idof in tabledata.dofmap]
            value = L.Sum(values)
            code = [
                L.VariableDecl("const double", access, value)
                ]
        elif inline:
            # Inlined version (we know this is bounded by a small number)
            dof_access = self.symbols.domain_dofs_access(gdim, num_scalar_dofs,
                                                         mt.restriction)
            # Inlined loop to accumulate linear combination of dofs and tables
            value = L.Sum([dof_access[idof] * FE[i]
                           for i, idof in enumerate(tabledata.dofmap)])
            code = [
                L.VariableDecl("const double", access, value)
                ]
        else:  # TODO: Make an option to test this version for performance
            # Assuming contiguous dofmap here
            assert len(tabledata.dofmap) == end - begin

            # Generated loop version:
            ic = self.symbols.coefficient_dof_sum_index()
            dof_access = self.symbols.domain_dof_access(ic + begin, mt.flat_component,
                                                        gdim, num_scalar_dofs,
                                                        mt.restriction)

            # Loop to accumulate linear combination of dofs and tables
            code = [
                L.VariableDecl("double", access, 0.0),
                L.ForRange(ic, 0, end - begin,
                           body=[L.AssignAdd(access, dof_access * FE[ic])])
                ]

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
        if self.integral_type in custom_integral_types:
            # FIXME: Jacobian may need adjustment for custom_integral_types
            if mt.local_derivatives:
                error("FIXME: Jacobian in custom integrals is not implemented.")
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
        expr = L.Conditional(L.EQ(co, L.LiteralInt(1)),
                             L.LiteralFloat(-1.0), L.LiteralFloat(+1.0))
        code = [
            L.VariableDecl("const double", access, expr)
            ]
        return code


    def _expect_table(self, e, mt, tabledata, num_points, access):
        "These quantities refer to constant tables defined in ufc_geometry.h."
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
        "These quantities refer to coordinate_dofs"
        # TODO: Generate more efficient inline code for Max/MinCell/FacetEdgeLength
        #       and CellDiameter here rather than lowering these quantities?
        return []
    cell_vertices = _expect_physical_coords
    cell_edge_vectors = _expect_physical_coords
    facet_edge_vectors = _expect_physical_coords


    def _expect_symbolic_lowering(self, e, mt, tabledata, num_points, access):
        "These quantities are expected to be replaced in symbolic preprocessing."
        error("Expecting {0} to be replaced in symbolic preprocessing.".format(type(e)))
    facet_normal = _expect_symbolic_lowering
    cell_normal = _expect_symbolic_lowering
    jacobian_inverse = _expect_symbolic_lowering
    jacobian_determinant = _expect_symbolic_lowering
    facet_jacobian = _expect_symbolic_lowering
    facet_jacobian_inverse = _expect_symbolic_lowering
    facet_jacobian_determinant = _expect_symbolic_lowering
