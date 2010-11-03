__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

from ufl.algorithms import preprocess

error_control_base = """
  class %(class_name)s: public dolfin::ErrorControl
  {

   public:

   %(class_name)s(const dolfin::Form& a, const dolfin::Form& L,
                  dolfin::GoalFunctional& M)
      : dolfin::%(class_name)s(a, L, M)

   {

     std::string name;

     // Create dual forms with function spaces from a
     _a_star.reset(new %(a_star)s(a.function_space(1), a.function_space(0)));
     _L_star.reset(new %(L_star)s(a.function_space(1)));

     %(attach_a_star)s
     %(attach_L_star)s

     const dolfin::FunctionSpace& V(*a.function_space(0));
     const dolfin::Mesh& mesh(V.mesh());

     // Create extrapolation space
     _E.reset(new %(E_space)s(mesh));

     // Create residual functional (for use with error estimate)
     _residual.reset(new %(residual)s(mesh));

     %(attach_residual)s

     // Create bilinear and linear form for computing cell residual R_T
     _DG_k.reset(new %(DGk_space)s(mesh));
     _a_R_T.reset(new %(a_R_T)s(*_DG_k, *_DG_k));
     _L_R_T.reset(new %(L_R_T)s(*_DG_k));

     // Initialize bubble function
     _B.reset(new %(Bubble_space)s(mesh));
     _b_T.reset(new dolfin::Function(*_B));
     _b_T->vector() = 1.0;

     // Attach bubble function to _a_R_T and _L_R_T
     _a_R_T->set_coefficient(0, *_b_T);
     _L_R_T->set_coefficient(%(L_R_T_b_T)d, *_b_T);

     %(attach_L_R_T)s

     // Create bilinear and linear form for computing facet residual R_dT
     _a_R_dT.reset(new %(a_R_dT)s(*_DG_k, *_DG_k));
     _L_R_dT.reset(new %(L_R_dT)s(*_DG_k));

     %(attach_L_R_dT)s

     // Initialize cone space
     _C.reset(new %(Cone_space)s(mesh));

     // Create error indicator form
     _DG_0.reset(new %(DG0_space)s(mesh));
     _eta_T.reset(new %(eta_T)s(*_DG_0));

   }

  ~%(class_name)s()
   {
     // Do nothing. Boost takes care of deletion
   }

  };

"""

typedefs = """
  typedef %(a)s BilinearForm;
  typedef %(L)s LinearForm;
  typedef %(M)s GoalFunctional;
  typedef %(a)s::TestSpace TestSpace;
  typedef %(a)s::TrialSpace TrialSpace;
"""

def generate_attach_snippet(to, from_form):

    snippet = """// Attach coefficients in %(to)s
    for (dolfin::uint i = 0; i < %(from_form)s.num_coefficients(); i++)
    {
      name = %(from_form)s.coefficient_name(i);
      std::cout << "Attaching coefficient named: " << name;
      std::cout << " to %(to)s" << std::endl;

      %(to)s->set_coefficient(name, %(from_form)s.coefficient(i));
    }""" % {"to": to, "from_form": from_form}

    return snippet

def generate_error_control_wrapper(prefix, maps):
    code = error_control_base % maps
    defs = typedefs % maps

    return (code, defs)

def generate_wrapper_maps(ec_forms, forms):
    """
    Form_0 : a_star
    Form_1 : L_star
    Form_2 : a_R_T
    Form_3 : L_R_T
    Form_4 : a_R_dT
    Form_5 : L_R_dT
    Form_6 : eta_h
    Form_7 : eta_T

    Form_8 : a
    Form_9 : L
    Form_10: M

    or

    Form_8 : F
    Form_9 : M

    Ez_h will be the last coefficient defined in Form_6

    """
    d = len(ec_forms) # Number of forms from the error control

    assert (d == 8), "Hm."
    (a_star, L_star, a_R_T, L_R_T, a_R_dT, L_R_dT, eta_h, eta_T) = ec_forms
    #(a_star, L_star, a_R_T, L_R_T) = ec_forms

    residual_snippet = generate_attach_snippet("_residual", "a") + "\n\t" + \
                       generate_attach_snippet("_residual", "L")
    L_R_T_snippet = generate_attach_snippet("_L_R_T", "a") + "\n\t" + \
                    generate_attach_snippet("_L_R_T", "L")
    L_R_dT_snippet = generate_attach_snippet("_L_R_dT", "a") + "\n\t" + \
                     generate_attach_snippet("_L_R_dT", "L")

    eta_h = preprocess(eta_h)
    print "eta_h = ", eta_h

    num_Ez_h = eta_h.form_data().num_coefficients - 1
    L_R_T = preprocess(L_R_T)
    num_b_T_in_L_R_T = L_R_T.form_data().num_coefficients - 1

    maps = {"class_name": "ErrorControl",
            "a_star": "Form_%d" % 0,
            "L_star": "Form_%d" % 1,
            "DGk_space": "Form_%d::TestSpace" % 3,
            "a_R_T": "Form_%d" % 2,
            "Bubble_space": "Form_%d_FunctionSpace_%d" % (2, 2),
            "L_R_T": "Form_%d" % 3,
            "a_R_dT": "Form_%d" % 4,
            "Cone_space": "Form_%d_FunctionSpace_%d" % (4, 2),
            "L_R_dT": "Form_%d" % 5,
            "residual": "Form_%d" % 6,
            "E_space": "Form_%d_FunctionSpace_%s" % (6, num_Ez_h),
            "eta_T": "Form_%d" % 7,
            "DG0_space": "Form_%d::TestSpace" % 7,
            "a": "Form_%d" % d,
            "L": "Form_%d" % (d+1),
            "M": "Form_%d" % (d+2),
            "L_R_T_b_T": num_b_T_in_L_R_T,
            "attach_a_star": generate_attach_snippet("_a_star", "a"),
            "attach_L_star": generate_attach_snippet("_L_star", "M"),
            "attach_residual": residual_snippet,
            "attach_L_R_T": L_R_T_snippet,
            "attach_L_R_dT": L_R_dT_snippet
            }

    return maps
def hack(lines, prefix):

    code = "".join(lines)

    n = 10

    (before, after) = code.split("class Form_%d:" % n)

    #print "after = ", after

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
