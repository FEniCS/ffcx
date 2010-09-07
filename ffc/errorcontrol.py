__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2010-09-07

from ufl.algorithms.analysis import extract_elements, extract_unique_elements
from ufl import adjoint, FiniteElement, MixedElement, Coefficient, action

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


def _check_input(forms):

    # At least check that we get three forms
    assert len(forms) == 3, "Not correct number of forms"

def _extract_forms(forms):

    # Extract separate forms (in a more robust way than this).
    (a, L, M) = forms

    return (a, L, M)

def create_dual_forms(a, L, M):

    # Optimal solution:
    # u = extract_trial_function(a)
    # a_star = adjoint(derivative(a - L, u))
    # L_star = derivative(M, u)

    # Fudged linear case for now
    a_star = adjoint(a)
    L_star = M

    return (a_star, L_star)

def create_extrapolation_space(L):
    """ The extrapolation space is a higher order version of the
    primal test space."""

    # Extract primal test space
    V = extract_unique_elements(L)[0]

    # Increase order and return
    return increase_order(V)

def generate_error_control_stuff(forms):

    # Check input
    _check_input(forms)

    # Extract forms
    (a, L, M) = _extract_forms(forms)

    # Create residual form (Probably need to name the coefficient
    # here, in order to be able to find it again.)
    residual = L - action(a)

    # Create forms for dual problem
    (a_star, L_star) = create_dual_forms(a, L, M)

    # Generate higher order element for extrapolation
    E = create_extrapolation_space(L)

    # Collect forms and elements to be compiled
    forms = (residual, a_star, L_star)
    elements = (E, )

    return (forms, elements)
