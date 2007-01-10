# Code generation format strings for UFC (Unified Form-assembly Code) v. 1.0.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenics.org/) 2006.

function_combined = """\
/// This class defines the interface for a general tensor-valued function.

class %(classname)s: public ufc::function
{%(members)s
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

  /// Evaluate function in the cell at the point
  /// x = (coordinates[0], coordinates[1], ...)
  void evaluate(double* values,
                const double* coordinates,
                const ufc::cell& c) const
  {
%(evaluate)s
  }

};
"""

function_header = """\
/// This class defines the interface for a general tensor-valued function.

class %(classname)s: public ufc::function
{%(members)s
public:

  /// Constructor
  %(classname)s();

  /// Destructor
  ~%(classname)s();

  /// Evaluate function in the cell at the point
  /// x = (coordinates[0], coordinates[1], ...)
  void evaluate(double* values,
                const double* coordinates,
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

/// Evaluate function in the cell at the point
/// x = (coordinates[0], coordinates[1], ...)
void %(classname)s::evaluate(double* values,
                             const double* coordinates,
                             const ufc::cell& c) const
{
%(evaluate)s
}
"""
