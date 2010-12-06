__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2010-12-06

from ufl import FiniteElement, MixedElement, Coefficient, TrialFunction, TestFunction
from ufl import adjoint, action, replace, inner, dx, ds, dS, avg, derivative
from ufl.algorithms.analysis import extract_elements, extract_unique_elements, extract_arguments

def change_regularity(element, family):
    """
    For a given function space, return the corresponding space with
    the finite elements specified by 'family'. Possible families
    are the families supported by the form compiler
    """

    # MR: This belongs in UFL
    n = element.num_sub_elements()
    if n > 0:
        subs = element.sub_elements()
        return MixedElement([change_regularity(subs[i], family)
                             for i in range(n)])
    shape = element.value_shape()
    if not shape:
        return FiniteElement(family, element.cell(), element.degree())

    return MixedElement([FiniteElement(family, element.cell(), element.degree())
                               for i in range(shape[0])])

def tear(V):
    "For a given space, return the corresponding discontinuous space."
    W = change_regularity(V, "DG")
    return W

def increase_order(element):
    "Return element of same family, but a polynomial degree higher."

    # MR: This belongs in UFL.
    n = element.num_sub_elements()
    if n > 0:
        subs = element.sub_elements()
        return MixedElement([increase_order(subs[i]) for i in range(n)])

    if element.family() == "Real":
        return element

    return FiniteElement(element.family(), element.cell(), element.degree()+1)

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
