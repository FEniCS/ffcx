# -*- coding: utf-8 -*-
"""
This module provides an abstract ErrorControlGenerator class for
generating forms required for goal-oriented error control and a
realization of this: UFLErrorControlGenerator for handling pure UFL
forms.
"""

# Copyright (C) 2010 Marie E. Rognes
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.

from ufl import inner, dx, ds, dS, avg, replace, action

__all__ = ["ErrorControlGenerator", "UFLErrorControlGenerator"]


class ErrorControlGenerator:

    def __init__(self, module, F, M, u):
        """
        *Arguments*

            module (Python module)
               The module to use for specific form manipulations
               (typically ufl or dolfin)

            F (tuple or Form)
               tuple of (bilinear, linear) forms or linear form

            M (Form)
               functional or linear form

            u (Coefficient)
              The coefficient considered as the unknown.
        """
        # Store module
        self.module = module

        # Store solution Coefficient/Function
        self.u = u

        # Extract the lhs (bilinear form), rhs (linear form), goal
        # (functional), weak residual (linear form)
        linear_case = (isinstance(F, (tuple, list)) and len(F) == 2)
        if linear_case:
            self.lhs, self.rhs = F
            try:
                self.goal = action(M, u)
            except Exception:
                self.goal = M    # Allow functionals as input as well
            self.weak_residual = self.rhs - action(self.lhs, u)
        else:
            self.lhs = self.module.derivative(F, u)
            self.rhs = F
            self.goal = M
            self.weak_residual = - F

        # At least check that final forms have correct rank
        assert(len(self.lhs.arguments()) == 2)
        assert(len(self.rhs.arguments()) == 1)
        assert(len(self.goal.arguments()) == 0)
        assert(len(self.weak_residual.arguments()) == 1)

        # Get the domain
        self.domain = self.weak_residual.ufl_domain()

        # Store map from identifiers to names for forms and generated
        # coefficients
        self.ec_names = {}

        # Use predefined names for the forms in the primal problem
        self.ec_names[id(self.lhs)] = "lhs"
        self.ec_names[id(self.rhs)] = "rhs"
        self.ec_names[id(self.goal)] = "goal"

        # Initialize other required data
        self.initialize_data()

    def initialize_data(self):
        """
        Initialize specific data
        """
        msg = """ErrorControlGenerator acts as an abstract
        class. Subclasses must overload the initialize_data() method
        and provide a certain set of variables. See
        UFLErrorControlGenerator for an example."""
        raise NotImplementedError(msg)

    def generate_all_error_control_forms(self):
        """
        Utility function for generating all (8) forms required for
        error control in addition to the primal forms
        """
        # Generate dual forms
        (a_star, L_star) = self.dual_forms()

        # Generate forms for computing strong cell residual
        (a_R_T, L_R_T) = self.cell_residual()

        # Generate forms for computing strong facet residuals
        (a_R_dT, L_R_dT) = self.facet_residual()

        # Generate form for computing error estimate
        eta_h = self.error_estimate()

        # Generate form for computing error indicators
        eta_T = self.error_indicators()

        # Paranoid checks added after introduction of multidomain features in ufl:
        for i, form in enumerate((a_star, L_star, eta_h, a_R_T, L_R_T, a_R_dT, L_R_dT, eta_T)):
            assert form.ufl_domain() is not None

        # Return all generated forms in CERTAIN order matching
        # constructor of dolfin/adaptivity/ErrorControl.h
        return (a_star, L_star, eta_h, a_R_T, L_R_T, a_R_dT, L_R_dT, eta_T)

    def primal_forms(self):
        """
        Return primal forms in order (bilinear, linear, functional)
        """
        return self.lhs, self.rhs, self.goal

    def dual_forms(self):
        """
        Generate and return (bilinear, linear) forms defining linear
        dual variational problem
        """
        a_star = self.module.adjoint(self.lhs)
        L_star = self.module.derivative(self.goal, self.u)
        return (a_star, L_star)

    def cell_residual(self):
        """
        Generate and return (bilinear, linear) forms defining linear
        variational problem for the strong cell residual
        """
        # Define trial and test functions for the cell residuals on
        # discontinuous version of primal trial space
        R_T = self.module.TrialFunction(self._dV)
        v = self.module.TestFunction(self._dV)

        # Extract original test function in the weak residual
        v_h = self.weak_residual.arguments()[0]

        # Define forms defining linear variational problem for cell
        # residual
        v_T = self._b_T * v
        a_R_T = inner(v_T, R_T) * dx(self.domain)
        L_R_T = replace(self.weak_residual, {v_h: v_T})

        return (a_R_T, L_R_T)

    def facet_residual(self):
        """
        Generate and return (bilinear, linear) forms defining linear
        variational problem for the strong facet residual(s)
        """
        # Define trial and test functions for the facet residuals on
        # discontinuous version of primal trial space
        R_e = self.module.TrialFunction(self._dV)
        v = self.module.TestFunction(self._dV)

        # Extract original test function in the weak residual
        v_h = self.weak_residual.arguments()[0]

        # Define forms defining linear variational problem for facet
        # residual
        v_e = self._b_e * v
        a_R_dT = ((inner(v_e('+'), R_e('+')) +
                   inner(v_e('-'), R_e('-'))) * dS(self.domain) +
                  inner(v_e, R_e) * ds(self.domain))
        L_R_dT = (replace(self.weak_residual, {v_h: v_e}) -
                  inner(v_e, self._R_T) * dx(self.domain))

        return (a_R_dT, L_R_dT)

    def error_estimate(self):
        """
        Generate and return functional defining error estimate
        """
        # Error estimate defined as r(Ez_h):
        eta_h = action(self.weak_residual, self._Ez_h)
        return eta_h

    def error_indicators(self):
        """
        Generate and return linear form defining error indicators
        """
        # Extract these to increase readability
        R_T = self._R_T
        R_dT = self._R_dT
        z = self._Ez_h
        z_h = self._z_h

        # Define linear form for computing error indicators
        v = self.module.TestFunction(self._DG0)
        eta_T = (v * inner(R_T, z - z_h) * dx(self.domain) +
                 avg(v)*(inner(R_dT('+'), (z - z_h)('+')) +
                         inner(R_dT('-'), (z - z_h)('-'))) * dS(self.domain) +
                 v * inner(R_dT, z - z_h) * ds(self.domain))

        return eta_T


class UFLErrorControlGenerator(ErrorControlGenerator):

    """
    This class provides a realization of ErrorControlGenerator for use
    with pure UFL forms
    """

    def __init__(self, F, M, u):
        """
        *Arguments*

            F (tuple or Form)
               tuple of (bilinear, linear) forms or linear form

            M (Form)
               functional or linear form

            u (Coefficient)
              The coefficient considered as the unknown.
        """

        ErrorControlGenerator.__init__(self, __import__("ufl"), F, M, u)

    def initialize_data(self):
        """
        Extract required objects for defining error control
        forms. This will be stored, reused and in particular named.
        """
        # Developer's note: The UFL-FFC-DOLFIN--PyDOLFIN toolchain for
        # error control is quite fine-tuned. In particular, the order
        # of coefficients in forms is (and almost must be) used for
        # their assignment. This means that the order in which these
        # coefficients are defined matters and should be considered
        # fixed.

        from ufl import FiniteElement, FunctionSpace, Coefficient
        from ufl.algorithms.elementtransformations import tear, increase_order

        # Primal trial element space
        self._V = self.u.ufl_function_space()

        # Extract domain
        domain = self.u.ufl_domain()

        # Primal test space == Dual trial space
        Vhat = self.weak_residual.arguments()[0].ufl_function_space()

        # Discontinuous version of primal trial element space
        self._dV = FunctionSpace(domain, tear(self._V.ufl_element()))

        # Extract geometric dimension
        gdim = domain.geometric_dimension()

        # Coefficient representing improved dual
        E = FunctionSpace(domain, increase_order(Vhat.ufl_element()))
        self._Ez_h = Coefficient(E)
        self.ec_names[id(self._Ez_h)] = "__improved_dual"

        # Coefficient representing cell bubble function
        Belm = FiniteElement("B", domain.ufl_cell(), gdim + 1)
        Bfs = FunctionSpace(domain, Belm)
        self._b_T = Coefficient(Bfs)
        self.ec_names[id(self._b_T)] = "__cell_bubble"

        # Coefficient representing strong cell residual
        self._R_T = Coefficient(self._dV)
        self.ec_names[id(self._R_T)] = "__cell_residual"

        # Coefficient representing cell cone function
        Celm = FiniteElement("DG", domain.ufl_cell(), gdim)
        Cfs = FunctionSpace(domain, Celm)
        self._b_e = Coefficient(Cfs)
        self.ec_names[id(self._b_e)] = "__cell_cone"

        # Coefficient representing strong facet residual
        self._R_dT = Coefficient(self._dV)
        self.ec_names[id(self._R_dT)] = "__facet_residual"

        # Define discrete dual on primal test space
        self._z_h = Coefficient(Vhat)
        self.ec_names[id(self._z_h)] = "__discrete_dual_solution"

        # Piecewise constants for assembling indicators
        self._DG0 = FunctionSpace(domain, FiniteElement("DG",
                                                        domain.ufl_cell(), 0))
