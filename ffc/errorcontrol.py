__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2010-10-13

from ufl.algorithms.analysis import extract_elements, extract_unique_elements, extract_arguments
from ufl import FiniteElement, MixedElement, Coefficient, TrialFunction, TestFunction
from ufl import adjoint, action, replace, inner, dx, ds, dS, avg

from ffc.log import info, error
from ffc.compiler import compile_form
from ffc.errorcontrolwrappers import generate_error_control_wrapper


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

def tear(element):
    """
    For a given element space, return the corresponding discontinuous
    space
    """
    W = change_regularity(element, "DG")
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

def create_cell_residual_forms(a, L, u_h):

    elements = extract_elements(a)

    # Define residual as linear form
    r = L - action(a, u_h)

    # Establish space for bubble
    cell = elements[0].cell()
    Bubble = FiniteElement("B", cell, cell.geometric_dimension()+1)

    # Define bubble
    b_T = Coefficient(Bubble)

    # Tear trial space
    DG = tear(elements[1])
    R_T = TrialFunction(DG)
    v = TestFunction(DG)

    # Define forms
    v_T = b_T*v
    a_R_T = inner(v_T, R_T)*dx
    L_R_T = replace(r, {extract_arguments(r)[0]: v_T})

    forms = (a_R_T, L_R_T)
    names = {"b_T": b_T}
    return (forms, names)


def create_facet_residual_forms(a, L, u):

    # Define weak residual (linear form)
    r = L - action(a, u)

    # Pick trial space of a as space of residual representation,
    elements = extract_elements(a)
    DG = tear(elements[1])

    # Establish cone function(s)
    cell = elements[0].cell()
    C = FiniteElement("DG", cell, cell.geometric_dimension())
    b_e = Coefficient(C)

    R_T = Coefficient(DG)
    R_e = TrialFunction(DG)
    v = TestFunction(DG)
    v_e = b_e*v

    a_R_dT = (inner(v_e('+'), R_e('+')) + inner(v_e('-'), R_e('-')))*dS \
             + inner(v_e, R_e)*ds
    L_R_dT = replace(r, {extract_arguments(r)[0]: v_e}) - inner(v_e, R_T)*dx

    forms = (a_R_dT, L_R_dT)
    names = {"b_e": b_e, "R_T": R_T}

    return (forms, names)

def create_error_indicator_form(a, L, z):

    elements = extract_elements(a)

    # R_T should live in a's trial space
    R_T = Coefficient(elements[1])

    # So should kind of R_dT
    R_dT = Coefficient(elements[1])

    # Interpolated dual extrapolation should live in dual trial space
    # aka primal test space
    z_h = Coefficient(elements[0])

    # Use test function on DG to localize
    DG = FiniteElement("DG", elements[1].cell(), 0)
    v = TestFunction(DG)

    # Define indicator form
    eta_T = v*inner(R_T, z - z_h)*dx \
            + avg(v)*(inner(R_dT('+'), (z - z_h)('+'))
                      + inner(R_dT('-'), (z - z_h)('-')))*dS \
            + v*inner(R_dT, z - z_h)*ds

    names = {"eta_T": eta_T, "R_T": R_T, "R_dT": R_dT, "z_h": z_h}

    return (eta_T, names)

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

    # Create bilinear and linear forms for cell residual
    (forms, R_T_names) = create_cell_residual_forms(a, L, u_h)
    (a_R_T, L_R_T) = forms

    # Create bilinear and linear forms for facet residual
    ((a_R_dT, L_R_dT), R_dT_names) = create_facet_residual_forms(a, L, u_h)

    # Create linear form for error indicators
    (eta_T, eta_T_names) = create_error_indicator_form(a, L, Ez_h)

    # Collect forms and elements to be compiled
    forms = (a_star, L_star, residual, a_R_T, L_R_T, a_R_dT, L_R_dT, eta_T)

    # Add names to object names
    names = {}
    names[id(a_star)] = "a_star"
    names[id(L_star)] = "L_star"
    names[id(residual)] = "residual"
    names[id(Ez_h)] = "Ez_h"
    names[id(u_h)] = "u_h"
    names[id(a_R_T)] = "a_R_T"
    names[id(L_R_T)] = "L_R_T"
    names[id(a_R_dT)] = "a_R_dT"
    names[id(L_R_dT)] = "L_R_dT"

    for (name, form) in eta_T_names.iteritems():
        names[id(form)] = name

    for (name, form) in R_T_names.iteritems():
        names[id(form)] = name

    for (name, form) in R_dT_names.iteritems():
        names[id(form)] = name

    return (forms, names)

def hack(lines, prefix):

    code = "".join(lines)

    n = 10

    (before, after) = code.split("class Form_%d:" % n)

    print "after = ", after

    after = after.replace(" public dolfin::Form",
                          "class Form_%d: public dolfin::GoalFunctional" % n)
    after = after.replace("dolfin::Form", "dolfin::GoalFunctional")

    updating = """
    virtual void update_ec(const dolfin::Form& a, const dolfin::Form& L)
    {
      // Update self
      ec.reset(new %s::ErrorControl(a, L, *this));
    }
    };
    """ % prefix

    after = after.replace("};", updating)

    return (before, after)

def write_code(prefix, ec_code, typedefs):

    # Append code to above file (must fix #endif)
    file = open(prefix + ".h", "r")
    lines = file.readlines()
    lines = lines[:-3]
    file.close()

    file = open(prefix + ".h", "w")
    # meg: Not proud moment. No need for yelling.
    (before, after) = hack(lines, prefix)
    file.write(before)
    file.write(ec_code)
    file.write(after)
    file.write(typedefs)
    file.write("}\n#endif\n");
    file.close()


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
    (ec_code, typedefs) = generate_error_control_wrapper(prefix)

    print "-"*80
    print "- Wrapper code - "
    print "-"*80
    print ec_code
    print "-"*80

    write_code(prefix, ec_code, typedefs)

    return 0

