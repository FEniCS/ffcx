__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2010-11-03

from ufl import Coefficient
from ufl.algorithms import preprocess
from ffc.log import info, error
from ffc.compiler import compile_form
from ffc.formmanipulations import *
from ffc.errorcontrolwrappers import *

def _check_input(forms, object_names):
    """
    Can handle three variants of forms:

    (F, M): F has rank 1 and M has rank 0
    (F, dF, M): F has rank 1, dF has rank 2 and M has rank 0
    (a, L, M): a has rank 2, L has rank a and M has rank 1

    """
    # FIXME: Add some more checks here...
    assert (len(forms) == 3 or len(forms) == 2), "Not correct number of forms"

def generate_error_control(forms, object_names):

    info("Generating additionals")

    _check_input(forms, object_names)

    ec_names = {}

    # Extract unknown (or None if not defined)
    unknown = object_names.get("unknown", None)

    # Generate dual forms
    a_star, L_star = generate_dual_forms(forms, unknown)

    # Extract trial space and generate extrapolation space from this
    V = extract_arguments(L_star)[0].element()
    E = increase_order(V)

    # If unknown is undefined, define discrete solution as coefficient
    # on trial element
    if unknown is None:
        unknown = Coefficient(V)
        ec_names[id(unknown)] = "the_discrete_solution"

    # Create coefficient for extrapolated dual
    Ez_h = Coefficient(E)

    # Create weak residual
    weak_residual = generate_weak_residual(forms[:-1], unknown)

    # Create cell residual forms
    a_R_T, L_R_T = generate_cell_residual(weak_residual)

    # Create facet residual forms
    a_R_dT, L_R_dT = generate_facet_residual(weak_residual)

    # Generate error estimate (residual) (# FIXME: Add option here)
    eta_h = action(weak_residual, Ez_h)

    print "eta_h = ", eta_h

    a = preprocess(eta_h)
    print "form_data = ", a.form_data()

    # Generate error indicators (# FIXME: Add option here)
    eta_T = generate_error_indicator(weak_residual, E, Q=None, Eh=None)

    ec_forms = (a_star, L_star, a_R_T, L_R_T, a_R_dT, L_R_dT, eta_h, eta_T)

    return (ec_forms, ec_names)

def compile_with_error_control(forms, object_names, prefix, parameters):

    ec_forms, ec_names = generate_error_control(forms, object_names)

    # Add object names generated for error control
    assert not (set(object_names.values()) & set(ec_names.values())), \
           "Same names used in code generation for error control. Trouble!"

    for (name, value) in ec_names.iteritems():
        object_names[name] = value

    print "\n\nObject names:"
    for (name, value) in object_names.iteritems():
        print (name, value)

    # Compile all forms
    all_forms = ec_forms + tuple(forms)
    compile_form(all_forms, object_names, prefix, parameters)

    maps = generate_wrapper_maps(ec_forms, forms)

    # Generate error_control DOLFIN wrapper
    (ec_code, typedefs) = generate_error_control_wrapper(prefix, maps)

    # Write code
    write_code(prefix, ec_code, typedefs)

    return 0



# def generate_linear_ec_forms(forms):

#     # Extract forms
#     #(a, L, M) = _extract_forms(forms)

#     # Create bilinear and linear forms for dual problem
#     (a_star, L_star) = create_dual_forms(forms)

#     # Define discrete solution as coefficient on trial element (NB!)
#     u_h = Coefficient(extract_elements(a)[1])

#     # Create finite element for extrapolation
#     E = create_extrapolation_space(L)

#     # Create coefficient for extrapolated dual
#     Ez_h = Coefficient(E)

#     # Create residual funcational
#     residual = action(L - action(a, u_h), Ez_h)

#     # Create bilinear and linear forms for cell residual
#     (forms, R_T_names) = create_cell_residual_forms(a, L, u_h)
#     (a_R_T, L_R_T) = forms

#     # Create bilinear and linear forms for facet residual
#     ((a_R_dT, L_R_dT), R_dT_names) = create_facet_residual_forms(a, L, u_h)

#     # Create linear form for error indicators
#     (eta_T, eta_T_names) = create_error_indicator_form(a, L, Ez_h)

#     # Collect forms and elements to be compiled
#     forms = (a_star, L_star, residual, a_R_T, L_R_T, a_R_dT, L_R_dT, eta_T)

#     # Add names to object names
#     names = {}
#     names[id(a_star)] = "a_star"
#     names[id(L_star)] = "L_star"
#     names[id(residual)] = "residual"
#     names[id(Ez_h)] = "Ez_h"
#     names[id(u_h)] = "u_h"
#     names[id(a_R_T)] = "a_R_T"
#     names[id(L_R_T)] = "L_R_T"
#     names[id(a_R_dT)] = "a_R_dT"
#     names[id(L_R_dT)] = "L_R_dT"

#     for (name, form) in eta_T_names.iteritems():
#         names[id(form)] = name

#     for (name, form) in R_T_names.iteritems():
#         names[id(form)] = name

#     for (name, form) in R_dT_names.iteritems():
#         names[id(form)] = name

#     return (forms, names)


