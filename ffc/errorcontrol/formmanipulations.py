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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC.  If not, see <http://www.gnu.org/licenses/>.
#
# Last changed: 2011-01-17

from ufl import adjoint, action, replace, inner, dx, ds, dS, avg, derivative
from ufl.algorithms.analysis import extract_arguments

__all__ = ["generate_dual_forms", "generate_weak_residual",
           "generate_cell_residual", "generate_facet_residual",
           "generate_error_indicator"]

def generate_dual_forms(forms, unknown, module):
    """
    Input:
    """
    (bilinear, linear, functional) = forms

    a_star = adjoint(bilinear)
    L_star = module.derivative(functional, unknown)

    return (a_star, L_star)

def generate_weak_residual(forms, u_h=None):

    if u_h is not None and len(forms) == 2:
        (a, L) = forms
        return L - action(a, u_h)

    return - forms

def generate_cell_residual(r, V_h, b_T, module):

    # Define trial and test functions
    R_T = module.TrialFunction(V_h)
    v = module.TestFunction(V_h)

    # Extract test argument from r for temporary use
    v_h = extract_arguments(r)[0]

    # Define forms
    v_T = b_T*v
    a_R_T = inner(v_T, R_T)*dx
    L_R_T = replace(r, {v_h: v_T})

    return (a_R_T, L_R_T)


def generate_facet_residual(r, V_h, b_e, R_T, module):

    v_h = extract_arguments(r)[0]

    # Define variational problem
    R_e = module.TrialFunction(V_h)
    v = module.TestFunction(V_h)

    # Define variational form
    v_e = b_e*v
    a_R_dT = (inner(v_e('+'), R_e('+')) + inner(v_e('-'), R_e('-')))*dS \
             + inner(v_e, R_e)*ds
    L_R_dT = replace(r, {v_h: v_e}) - inner(v_e, R_T)*dx

    return (a_R_dT, L_R_dT)

def generate_error_indicator(r, R_T, R_dT, z, z_h, v):

    # Define indicator form
    eta_T = v*inner(R_T, z - z_h)*dx \
            + avg(v)*(inner(R_dT('+'), (z - z_h)('+'))
                      + inner(R_dT('-'), (z - z_h)('-')))*dS \
            + v*inner(R_dT, z - z_h)*ds

    return eta_T
