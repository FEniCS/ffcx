__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"


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
     _a_R_T->set_coefficient("b_T", *_b_T);
     _L_R_T->set_coefficient("b_T", *_b_T);

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

def generate_error_control_wrapper(prefix):

    N = 0 # Assuming that forms are generated in order
    d = 8 # Number of forms from the error control

    residual_snippet = generate_attach_snippet("_residual", "a") + "\n\t" + \
                       generate_attach_snippet("_residual", "L")

    L_R_T_snippet = generate_attach_snippet("_L_R_T", "a") + "\n\t" + \
                    generate_attach_snippet("_L_R_T", "L")

    L_R_dT_snippet = generate_attach_snippet("_L_R_dT", "a") + "\n\t" + \
                     generate_attach_snippet("_L_R_dT", "L")

    maps = {"class_name": "ErrorControl",
            "a_star": "Form_%d" % (N),
            "L_star": "Form_%d" % (N+1),
            "E_space": "CoefficientSpace_%s" % "Ez_h",
            "residual": "Form_%d" % (N+2),
            "DGk_space": "Form_%d::TestSpace" % (N+3),
            "a_R_T": "Form_%d" % (N+3),
            "L_R_T": "Form_%d" % (N+4),
            "Bubble_space": "Form_%d::CoefficientSpace_%s" % (N+3, "b_T"),
            "Cone_space": "Form_%d::CoefficientSpace_%s" % (N+5, "b_e"),
            "a_R_dT": "Form_%d" % (N+5),
            "L_R_dT": "Form_%d" % (N+6),
            "DG0_space": "Form_%d::TestSpace" % (N+7),
            "eta_T": "Form_%d" % (N+7),
            "a": "Form_%d" %(N + d),
            "L": "Form_%d" %(N + d+1),
            "M": "Form_%d" %(N + d+2),
            "attach_a_star": generate_attach_snippet("_a_star", "a"),
            "attach_L_star": generate_attach_snippet("_L_star", "M"),
            "attach_residual": residual_snippet,
            "attach_L_R_T": L_R_T_snippet,
            "attach_L_R_dT": L_R_dT_snippet
            }

    code = error_control_base % maps
    defs = typedefs % maps

    return (code, defs)
