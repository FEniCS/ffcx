__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU LGPL version 3 or any later version"

base = """
class %(class_name)s: public dolfin::ErrorControl
{

 public:

  %(class_name)s(dolfin::Form& a, dolfin::Form& L, dolfin::Form& M)
    {

      // Create dual forms with function spaces from a
      _a_star = new %(a_star)s(a.function_space(1), a.function_space(0));
      _L_star = new %(L_star)s(a.function_space(1));

      // Create extrapolation space
      _E = new %(E_space)s(mesh);

      // Create residual functional (for use with error estimate)
      _residual = new %(residual)s(mesh);

    };

  ~%(class_name)s() {

    // delete stuff;

  };

};
"""

def generate_error_control_wrapper(prefix):

    N = 0 # Assuming that forms are generated in order
    code = base % {"class_name": prefix + "_ErrorControl",
                   "a_star": "Form_%d" % (N),
                   "L_star": "Form_%d" % (N+1),
                   "E_space": "CoefficientSpace_%s" % "Ez_h",
                   "residual": "Form_%d" % (N+2)}
    return code
