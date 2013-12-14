%module ufc

%{
#include "ufc.h"
// If using std::tr1::shared_ptr comment out this line
#include <boost/shared_ptr.hpp>
// Un comment this line to use std::tr1, only works with swig version  >= 1.3.37
//#include <tr1/memory>
%}

// Un comment these lines to use std::tr1, only works with swig version  >= 1.3.37
//#define SWIG_SHARED_PTR_NAMESPACE std
//#define SWIG_SHARED_PTR_SUBNAMESPACE tr1
%include <boost_shared_ptr.i>

// Ignore interface to ufc::cell that will not be available for the user
%ignore ufc::cell::entity_indices;

// Declare which classes should be stored using shared_ptr
%shared_ptr(ufc::cell_integral)
%shared_ptr(ufc::dofmap)
%shared_ptr(ufc::finite_element)
%shared_ptr(ufc::function)
%shared_ptr(ufc::form)
%shared_ptr(ufc::exterior_facet_integral)
%shared_ptr(ufc::interior_facet_integral)
%shared_ptr(ufc::point_integral)

%include <exception.i>

//-----------------------------------------------------------------------------
// Home brewed versions of the SWIG provided SWIG_AsVal(Type).
//-----------------------------------------------------------------------------
%fragment("Py_convert_uint", "header") {
  // A check for int and converter to uint
  SWIGINTERNINLINE bool Py_convert_uint(PyObject* in, std::size_t& value)
  {
    if (!(PyInt_Check(in) && PyInt_AS_LONG(in)>=0))
      return false;
    value = static_cast<std::size_t>(PyInt_AS_LONG(in));
    return true;
  }
}

//-----------------------------------------------------------------------------
// Out typemap (std::size_t)
//-----------------------------------------------------------------------------
%typemap(out) std::size_t
{
  // Typemap std::size_t
  $result = PyInt_FromLong(static_cast< long >($1));
}

//-----------------------------------------------------------------------------
// Typecheck and in typemap (std::size_t)
//-----------------------------------------------------------------------------
%typecheck(SWIG_TYPECHECK_INTEGER) std::size_t
{
  $1 = PyInt_Check($input) ? 1 : 0;
}

%typemap(in, fragment="Py_convert_uint") std::size_t
{
  if (!Py_convert_uint($input, $1))
    SWIG_exception(SWIG_TypeError, "expected positive 'int' for argument $argnum");
}



//-----------------------------------------------------------------------------
// Include the main header file
//-----------------------------------------------------------------------------
%include "ufc.h"

// Include code to generate a __swigversion__ attribute to the cpp module
// Add prefix to avoid naming problems with other modules
%inline %{
int ufc_swigversion() { return SWIGVERSION; }
%}

%pythoncode %{
__version__ = UFC_VERSION
del UFC_VERSION, UFC_VERSION_MAJOR, UFC_VERSION_MINOR

"""Code for adding swig version to ufc extension module."""
tmp = hex(ufc_swigversion())
__swigversion__ = "%d.%d.%d"%(tuple(map(int, [tmp[-5], tmp[-3], tmp[-2:]])))

del tmp, ufc_swigversion
%}
