
function_combined = """\
/// This class defines the interface for a general tensor-valued function.

class %(classname)s: public ufc::function
{
%(members)s
public:

  /// Constructor
  %(classname)s()
  {
%(constructor)s
  }

  /// Destructor
  ~%(classname)s()
  {
%(destructor)s
  }

  /// Evaluate the function at the point x = (x[0], x[1], ...) in the cell
  void evaluate(double* values,
                const double* x,
                const ufc::cell& c) const
  {
%(evaluate)s
  }

};
"""


function_header = """\
/// This class defines the interface for a general tensor-valued function.

class %(classname)s: public ufc::function
{
%(members)s
public:

  /// Constructor
  %(classname)s();

  /// Destructor
  ~%(classname)s();

  /// Evaluate the function at the point x = (x[0], x[1], ...) in the cell
  void evaluate(double* values,
                const double* x,
                const ufc::cell& c) const;

};
"""


function_implementation = """\
/// Constructor
%(classname)s::%(classname)s()
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Evaluate the function at the point x = (x[0], x[1], ...) in the cell
void %(classname)s::evaluate(double* values,
                             const double* x,
                             const ufc::cell& c) const
{
%(evaluate)s
}
"""

