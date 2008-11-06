"Templates for generating DOLFIN wrappers"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2008-11-06 -- 2008-11-06"
__copyright__ = "Copyright (C) 2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

dolfin_includes = """\
// DOLFIN wrappers

#include <dolfin/fem/Form.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/function/Coefficient.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>

"""

# Code for adding function space (reference version)
add_function_space_r = """\
    std::tr1::shared_ptr<const dolfin::FunctionSpace> _V%d(&V%d, dolfin::NoDeleter<const dolfin::FunctionSpace>());
    _function_spaces.push_back(_V%d);"""

# Code for adding function space (shared pointer version)
add_function_space_s = """\
    _function_spaces.push_back(V%d);"""

# Code for adding coefficient (reference version)
add_coefficient_r =  """\
    _coefficients.push_back(std::tr1::shared_ptr<const dolfin::Function>(static_cast<const dolfin::Function*>(0)));"""

# Code for adding coefficient (shared pointer version)
add_coefficient_s =  """\
    _coefficients.push_back(std::tr1::shared_ptr<const dolfin::Function>(static_cast<const dolfin::Function*>(0)));"""

# Code for coefficient class
coefficient_class = """\
class %s : public dolfin::Coefficient
{
public:

  // Constructor
  %s(dolfin::Form& form) : dolfin::Coefficient(form) {}

  // Destructor  
  ~%s() {}

  // Attach function to coefficient
  const %s& operator= (dolfin::Function& v)
  {
    attach(v);
    return *this;
  }

  /// Create function space for coefficient
  const dolfin::FunctionSpace* create_function_space() const
  {
    return new %sCoefficientSpace%d(form.mesh());
  }
  
  /// Return coefficient number
  dolfin::uint number() const
  {
    return %d;
  }
  
  /// Return coefficient name
  virtual std::string name() const
  {
    return "%s";
  }
  
};
"""

# Code for form class (including function spaces and coefficients)
form_class_vc = """\
class %s : public dolfin::Form
{
public:

  // Create form on given function space(s)
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Create form on given function space(s) (shared data)
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Create form on given function space(s) with given coefficient(s)
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Create form on given function space(s) with given coefficient(s) (shared data)
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Destructor
  ~%s() {}
%s
};

"""

# Code for form class (including only function spaces)
form_class_v = """\
class %s : public dolfin::Form
{
public:

  // Create form on given function space(s)
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Create form on given function space(s) (shared data)
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Destructor
  ~%s() {}
%s
};

"""

# Code for form class (including only functions)
form_class_c = """\
class %s : public dolfin::Form
{
public:

  // Create form
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Create form with given coefficient(s)
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Destructor
  ~%s() {}
%s
};

"""

# Code for form class (with no arguments)
form_class = """\
class %s : public dolfin::Form
{
public:

  // Create form
  %s(%s) : dolfin::Form()%s
  {
%s
  }

  // Destructor
  ~%s() {}
%s
};

"""
