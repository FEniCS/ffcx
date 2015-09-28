# Code generation format strings for UFC (Unified Form-assembly Code) v. 1.7.0dev.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2015

function_combined = """\
/// This class defines the interface for a general tensor-valued function.

class %(classname)s: public ufc::function
{%(members)s
public:

  /// Constructor
  %(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::function()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Evaluate function at given point in cell
  void evaluate(double * values,
                const double * coordinates,
                const ufc::cell& c) const final override
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
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  ~%(classname)s() override;

  /// Evaluate function at given point in cell
  void evaluate(double * values,
                const double * coordinates,
                const ufc::cell& c) const final override;

};
"""

function_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::function()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Evaluate function at given point in cell
void %(classname)s::evaluate(double * values,
                             const double * coordinates,
                             const ufc::cell& c) const
{
%(evaluate)s
}
"""
