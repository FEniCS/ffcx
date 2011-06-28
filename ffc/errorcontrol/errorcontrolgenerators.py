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
#
# Last changed: 2011-06-28

from ufl import inner, dx, ds, dS, avg, adjoint, replace, action
from ufl.algorithms.analysis import extract_arguments

__all__ = ["ErrorControlGenerator", "UFLErrorControlGenerator"]

class ErrorControlGenerator:

    def __init__(self, module, F, M, u):

        # Store module
        self.module = module

        # Store solution Coefficient/Function
        self.u = u

        # Extract the lhs (bilinear form), rhs (linear form), goal
        # (functional), weak residual (linear form)
        # FIXME: MER: Error checking is not for whimps.
        if (isinstance(F, (tuple, list)) and len(F) == 2):
            self.lhs, self.rhs = F
            try:
                self.goal = action(M, u)
            except:
                self.goal = M
            self.weak_residual = self.rhs - action(self.lhs, u)
        else:
            self.lhs = self.module.derivative(F, u)
            self.rhs = F
            self.goal = M
            self.weak_residual = - F

    def initialize_data(self):
        # FIXME: Improve robustness here.
        msg = "Please overload initialize_data and know what you are doing"
        raise NotImplementedError, msg

    def generate_all_error_control_forms(self):

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

        # Return all generated forms in sensible order
        return (a_star, L_star, eta_h, a_R_T, L_R_T, a_R_dT, L_R_dT, eta_T)

    def primal_forms(self):
        return self.lhs, self.rhs, self.goal

    def dual_forms(self):
        a_star = adjoint(self.lhs)
        L_star = self.module.derivative(self.goal, self.u)
        return (a_star, L_star)

    def cell_residual(self):

        # Define trial and test functions for the cell residuals on
        # discontinuous version of primal trial space
        R_T = self.module.TrialFunction(self._dV)
        v = self.module.TestFunction(self._dV)

        # Extract original test function in the weak residual
        v_h = extract_arguments(self.weak_residual)[0]

        # Define forms defining linear variational problem for cell
        # residual
        v_T = self._b_T*v
        a_R_T = inner(v_T, R_T)*dx
        L_R_T = replace(self.weak_residual, {v_h: v_T})

        return (a_R_T, L_R_T)

    def facet_residual(self):

        # Define trial and test functions for the facet residuals on
        # discontinuous version of primal trial space
        R_e = self.module.TrialFunction(self._dV)
        v = self.module.TestFunction(self._dV)

        # Extract original test function in the weak residual
        v_h = extract_arguments(self.weak_residual)[0]

        # Define forms defining linear variational problem for facet
        # residual
        v_e = self._b_e*v
        a_R_dT = ((inner(v_e('+'), R_e('+')) + inner(v_e('-'), R_e('-')))*dS
                  + inner(v_e, R_e)*ds)
        L_R_dT = (replace(self.weak_residual, {v_h: v_e})
                  - inner(v_e, self._R_T)*dx)

        return (a_R_dT, L_R_dT)

    def error_estimate(self):
        # eta_h = r(Ez_h)
        return action(self.weak_residual, self._Ez_h)

    def error_indicators(self):

        # Extract these to increase readability
        R_T = self._R_T
        R_dT = self._R_dT
        z = self._Ez_h
        z_h = self._z_h

        v = self.module.TestFunction(self._DG0)
        eta_T = (v*inner(R_T, z - z_h)*dx
                 + avg(v)*(inner(R_dT('+'), (z - z_h)('+'))
                           + inner(R_dT('-'), (z - z_h)('-')))*dS
                 + v*inner(R_dT, z - z_h)*ds)
        return eta_T

class UFLErrorControlGenerator(ErrorControlGenerator):
    def __init__(self, F, M, u):

        ErrorControlGenerator.__init__(self, __import__("ufl"), F, M, u)
        self.ec_names = {}
        self.initialize_data()

    def initialize_data(self):
        # Extract and store objects that will be shared between forms,
        # and in particular possibly be named

        from ufl import FiniteElement, Coefficient
        from ufl.algorithms.elementtransformations import tear, increase_order

        # Primal trial element space
        self._V = self.u.element()

        # Primal test space == Dual trial space
        Vhat = extract_arguments(self.weak_residual)[0].element()

        # Discontinuous version of primal trial element space
        self._dV = tear(self._V)

        # Extract cell and geometric dimension
        cell = self._V.cell()
        g_dim = cell.geometric_dimension()

        # Coefficient representing improved dual
        E = increase_order(Vhat)
        self._Ez_h = Coefficient(E)
        self.ec_names[id(self._Ez_h)] = "__improved_dual"

        # Coefficient representing cell bubble function
        B = FiniteElement("B", cell, g_dim + 1)
        self._b_T = Coefficient(B)
        self.ec_names[id(self._b_T)] = "__cell_bubble"

        # Coefficient representing strong cell residual
        self._R_T = Coefficient(self._dV)
        self.ec_names[id(self._R_T)] = "__cell_residual"

        # Coefficient representing cell cone function
        C = FiniteElement("DG", cell, g_dim)
        self._b_e = Coefficient(C)
        self.ec_names[id(self._b_e)] = "__cell_cone"

        # Coefficient representing strong facet residual
        self._R_dT = Coefficient(self._dV)
        self.ec_names[id(self._R_dT)] = "__facet_residual"

        # Define discrete dual on primal test space
        self._z_h = Coefficient(Vhat)
        self.ec_names[id(self._z_h)] = "__discrete_dual_solution"


        self._DG0 = FiniteElement("DG", cell, 0)

