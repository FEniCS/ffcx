__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

# Last changed: 2010-12-06

__all__ = ["generate_ec_generator_code", "generate_ec_typedefs",
           "write_code"]

ec_typedefs = """
typedef %(a)s BilinearForm;
typedef %(L)s LinearForm;
typedef %(M)s GoalFunctional;
typedef %(a)s::TestSpace TestSpace;
typedef %(a)s::TrialSpace TrialSpace;
"""

attach_snippet_base = """
    for (dolfin::uint i = 0; i < %(from)s.num_coefficients(); i++)
    {
      name = %(from)s.coefficient_name(i);
      // Don't attach discrete primal solution here (not computed.)
      if (name == "__discrete_primal_solution")
        continue;

      try {
        const uint coefficient_number = %(to)s->coefficient_number(name);
      } catch (...) {
        std::cout << "Attaching coefficient named: " << name << " to %(to)s";
        std::cout << " failed! But this might be expected." << std::endl;
        continue;
      }

      std::cout << "Attaching coefficient named: " << name;
      std::cout << " to %(to)s" << std::endl;
      %(to)s->set_coefficient(name, %(from)s.coefficient(i));
    }
    """

update_ec_base = """
  /// Initialize all error control forms, attach coefficients and
  /// (re-)set error control
  virtual void update_ec(const dolfin::Form& a, const dolfin::Form& L)
  {
    // This stuff is created here and shipped elsewhere
    boost::shared_ptr<dolfin::Form> a_star;           // Dual lhs
    boost::shared_ptr<dolfin::Form> L_star;           // Dual rhs
    boost::shared_ptr<dolfin::FunctionSpace> V_Ez_h;  // Extrapolation space
    boost::shared_ptr<dolfin::Function> Ez_h;         // Extrapolated dual
    boost::shared_ptr<dolfin::Form> residual;         // Residual (as functional)
    boost::shared_ptr<dolfin::FunctionSpace> V_R_T;   // Trial space for cell residual
    boost::shared_ptr<dolfin::Form> a_R_T;            // Cell residual lhs
    boost::shared_ptr<dolfin::Form> L_R_T;            // Cell residual rhs
    boost::shared_ptr<dolfin::FunctionSpace> V_b_T;   // Function space for cell bubble
    boost::shared_ptr<dolfin::Function> b_T;          // Cell bubble
    boost::shared_ptr<dolfin::FunctionSpace> V_R_dT;  // Trial space for facet residual
    boost::shared_ptr<dolfin::Form> a_R_dT;           // Facet residual lhs
    boost::shared_ptr<dolfin::Form> L_R_dT;           // Facet residual rhs
    boost::shared_ptr<dolfin::FunctionSpace> V_b_e;   // Function space for cell cone
    boost::shared_ptr<dolfin::Function> b_e;          // Cell cone
    boost::shared_ptr<dolfin::FunctionSpace> V_eta_T; // Function space for indicators
    boost::shared_ptr<dolfin::Form> eta_T;            // Indicator form

    // Some handy views
    const dolfin::FunctionSpace& Vhat(*(a.function_space(0))); // Primal test
    const dolfin::FunctionSpace& V(*(a.function_space(1)));    // Primal trial
    const dolfin::Mesh& mesh(V.mesh());
    std::string name;

    // Initialize dual forms
    a_star.reset(new %(a_star)s(V, Vhat));
    L_star.reset(new %(L_star)s(V));

    // Attach coefficients from a to a_star and from M to L_star
    %(attach_a_star)s
    %(attach_L_star)s

    // Initialize residual
    residual.reset(new %(residual)s(mesh));

    // Attach coefficients (from a and L) in residual
    %(attach_residual)s

    // Initialize extrapolation space and (fake) extrapolation
    V_Ez_h.reset(new %(V_Ez_h)s(mesh));
    Ez_h.reset(new dolfin::Function(V_Ez_h));
    std::cout << "Attaching (fake) __improved dual" << std::endl;
    residual->set_coefficient("__improved_dual", Ez_h);

    // Create bilinear and linear form for computing cell residual R_T
    V_R_T.reset(new %(V_R_T)s(mesh));
    a_R_T.reset(new %(a_R_T)s(V_R_T, V_R_T));
    L_R_T.reset(new %(L_R_T)s(V_R_T));

    // Initialize bubble and attach to a_R_T and L_R_T
    V_b_T.reset(new %(V_b_T)s(mesh));
    b_T.reset(new dolfin::Function(V_b_T));
    b_T->vector() = 1.0;

    // Attach coefficients (from a and L) to L_R_T
    %(attach_L_R_T)s

    // Attach bubble function to _a_R_T and _L_R_T
    std::cout << "Attaching __cell_bubble" << std::endl;
    a_R_T->set_coefficient("__cell_bubble", b_T);
    L_R_T->set_coefficient("__cell_bubble", b_T);

    // Create bilinear and linear form for computing facet residual R_dT
    V_R_dT.reset(new %(V_R_dT)s(mesh));
    a_R_dT.reset(new %(a_R_dT)s(V_R_dT, V_R_dT));
    L_R_dT.reset(new %(L_R_dT)s(V_R_dT));

    // Attach coefficients (from a and L) to L_R_dT
    %(attach_L_R_dT)s

    // Initialize (fake) cone and attach to a_R_dT and L_R_dT
    V_b_e.reset(new %(V_b_e)s(mesh));
    b_e.reset(new dolfin::Function(V_b_e));
    std::cout << "Attaching __cell_cone" << std::endl;
    a_R_dT->set_coefficient("__cell_cone", b_e);
    L_R_dT->set_coefficient("__cell_cone", b_e);

    // Create error indicator form
    V_eta_T.reset(new %(V_eta_T)s(mesh));
    eta_T.reset(new %(eta_T)s(V_eta_T));

    // Update error control
    _ec.reset(new dolfin::ErrorControl(a_star, L_star, residual,
                                       a_R_T, L_R_T, a_R_dT, L_R_dT, eta_T));

  }

"""

def generate_wrapper_maps():

    def attach(tos, froms):
        if isinstance(froms, tuple):
            return "\n".join([attach_snippet_base % {"to": to, "from": fro}
                              for (to, fro) in zip(tos, froms)])
        return attach_snippet_base % {"to": tos, "from": froms}

    maps = {"a_star":           "Form_%d" % 0,
            "L_star":           "Form_%d" % 1,
            "a_R_T":            "Form_%d" % 2,
            "L_R_T":            "Form_%d" % 3,
            "a_R_dT":           "Form_%d" % 4,
            "L_R_dT":           "Form_%d" % 5,
            "residual":         "Form_%d" % 6,
            "eta_T":            "Form_%d" % 7,
            "V_Ez_h":           "CoefficientSpace_%s" % "__improved_dual",
            "V_R_T":            "Form_%d::TestSpace" % 3,
            "V_b_T":            "CoefficientSpace_%s" % "__cell_bubble",
            "V_R_dT":           "Form_%d::TestSpace" % 5,
            "V_b_e":            "CoefficientSpace_%s" % "__cell_cone",
            "V_eta_T":          "Form_%d::TestSpace" % 7,
            "attach_a_star":    attach("a_star", "a"),
            "attach_L_star":    attach("L_star", "(*this)"),
            "attach_residual":  attach(("residual",)*2, ("a", "L")),
            "attach_L_R_T":     attach(("L_R_T",)*2, ("a", "L")),
            "attach_L_R_dT":    attach(("L_R_dT",)*2, ("a", "L"))
            }

    return maps

def generate_ec_generator_code():

    maps = generate_wrapper_maps()
    code = update_ec_base % maps
    return code

def generate_ec_typedefs(offset):

    maps = {"a": "Form_%d" % offset,
            "L": "Form_%d" % (offset+1),
            "M": "Form_%d" % (offset+2),
            }
    return ec_typedefs % maps

def write_code(prefix, class_name, extra_method, typedefs):

    # meg: Not proud moment. No need for yelling.

    # Append code to above file (must fix #endif)
    file = open(prefix + ".h", "r")
    lines = file.readlines()
    file.close()

    K = 5
    last = lines[-K:]
    lines = lines[:-K]

    code = "".join(lines)
    n = 10
    (before, after) = code.split("class Form_%d:" % 10)
    after = after.replace(" public dolfin::Form",
                          "class Form_%d: public %s" % (n, class_name))
    after = after.replace("dolfin::Form", class_name)
    after += extra_method

    file = open(prefix + ".h", "w")
    file.write(before)
    file.write(after)
    file.write("".join(last[0:2]))
    file.write(typedefs)
    file.write("".join(last[2:]))
    file.close()


