__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2010-09-08

from ufl.algorithms.analysis import extract_elements, extract_unique_elements, extract_arguments
from ufl import FiniteElement, MixedElement, Coefficient
from ufl import adjoint, action, replace

from ffc.log import info, error
from ffc.compiler import compile_form
from ffc.errorcontrolwrappers import generate_error_control_wrapper


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
    (u, v) = extract_arguments(a)
    L_star = replace(M, {u:v})

    return (a_star, L_star)

def create_extrapolation_space(L):
    """ The extrapolation space is a higher order version of the
    primal test space."""

    # Extract primal test space
    V = extract_unique_elements(L)[0]

    # Increase order and return
    return increase_order(V)

def generate_error_control_forms(forms):

    # Check input
    _check_input(forms)

    # Extract forms
    (a, L, M) = _extract_forms(forms)

    # Create bilinear and linear forms for dual problem
    (a_star, L_star) = create_dual_forms(a, L, M)

    # Define discrete solution as coefficient on trial element (NB!)
    u_h = Coefficient(extract_elements(a)[1])

    # Create finite element for extrapolation
    E = create_extrapolation_space(L)

    # Create coefficient for extrapolated dual
    Ez_h = Coefficient(E)

    # Create residual funcational
    residual = action(L - action(a, u_h), Ez_h)

    # Collect forms and elements to be compiled
    forms = (a_star, L_star, residual)

    # Add names to object names
    names = {}
    names[id(a_star)] = "a_star"
    names[id(L_star)] = "L_star"
    names[id(residual)] = "residual"
    names[id(Ez_h)] = "Ez_h"
    names[id(u_h)] = "u_h"

    return (forms, names)

def compile_with_error_control(forms, object_names, prefix, parameters):

    info("Generating additionals")
    (foos, names) = generate_error_control_forms(forms)

    # Check whether we use same names...
    if bool(set(names.values()) & set(object_names.values())):
        error("Same name used ... this can cause trouble")

    # Note: Not quite sure what to use this for yet.
    all_names = {}
    for k in names:
        all_names[k] = names[k]
    for k in object_names:
        all_names[k] = object_names[k]

    # Compile all forms
    all_forms = foos + tuple(forms)
    compile_form(all_forms, all_names, prefix, parameters)

    # Generate error_control DOLFIN wrapper
    code = generate_error_control_wrapper(prefix)

    print "-"*80
    print "- Wrapper code - "
    print "-"*80
    print code
    print "-"*80

    # Append code to above file (must fix #endif)
    file = open(prefix + ".h", "a")
    file.write(code)
    file.close()

    return 0

