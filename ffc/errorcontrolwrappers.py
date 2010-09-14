__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

base = """
  class %(class_name)s: public dolfin::ErrorControl
  {

   public:

   %(class_name)s(const dolfin::Form& a, const dolfin::Form& L, dolfin::Form& M)
      : dolfin::%(class_name)s(a, L, M)

   {

     std::string name;

     // Create dual forms with function spaces from a
     _a_star = new %(a_star)s(a.function_space(1), a.function_space(0));
     _L_star = new %(L_star)s(a.function_space(1));

     %(attach_a_star)s
     %(attach_L_star)s

     const dolfin::FunctionSpace& V(*a.function_space(0));
     const dolfin::Mesh& mesh(V.mesh());

     // Create extrapolation space
     _E = new %(E_space)s(mesh);

     // Create residual functional (for use with error estimate)
     _residual = new %(residual)s(mesh);

     %(attach_residual)s

     // Create bilinear and linear form for computing cell residual R_T
     _DG_k = new %(DGk_space)s(mesh);
     _a_R_T = new %(a_R_T)s(*_DG_k, *_DG_k);
     _L_R_T = new %(L_R_T)s(*_DG_k);

     // Create error indicator form
     _DG_0 = new %(DG0_space)s(mesh);
     _eta_T = new %(eta_T)s(*_DG_0);

   }

  ~%(class_name)s()
   {
     // delete stuff;
   }

  };

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
    d = 6 # Number of forms from the error control

    residual_snippet = generate_attach_snippet("_residual", "a") + "\n" + \
                       generate_attach_snippet("_residual", "L")

    code = base % {"class_name": "ErrorControl",
                   "a_star": "Form_%d" % (N),
                   "L_star": "Form_%d" % (N+1),
                   "E_space": "CoefficientSpace_%s" % "Ez_h",
                   "residual": "Form_%d" % (N+2),
                   "DGk_space": "Form_%d::TestSpace" % (N+3),
                   "a_R_T": "Form_%d" % (N+3),
                   "L_R_T": "Form_%d" % (N+4),
                   "DG0_space": "Form_%d::TestSpace" % (N+5),
                   "eta_T": "Form_%d" % (N+5),
                   "a": "Form_%d" %(N + d),
                   "L": "Form_%d" %(N + d+1),
                   "M": "Form_%d" %(N + d+2),
                   "attach_a_star": generate_attach_snippet("_a_star", "a"),
                   "attach_L_star": generate_attach_snippet("_L_star", "M"),
                   "attach_residual": residual_snippet}

    return code
