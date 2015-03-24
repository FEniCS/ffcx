%module ufc

%{
#include "ufc.h"
#include <memory>
%}

// Use std::shared_ptr
%include <std_shared_ptr.i>

// Ignore interface to ufc::cell that will not be available for the user
%ignore ufc::cell::entity_indices;

// Declare which classes should be stored using shared_ptr
%include "ufc_shared_ptr_classes.i"
%include <exception.i>

//-----------------------------------------------------------------------------
// Home brewed versions of the SWIG provided SWIG_AsVal(Type).
//-----------------------------------------------------------------------------
%fragment("Py_convert_uint", "header") {
  // A check for int and converter to uint
  SWIGINTERNINLINE bool Py_convert_uint(PyObject* in, std::size_t& value)
  {
%#if PY_MAJOR_VERSION >= 3
    if (!(PyLong_Check(in) && PyLong_AS_LONG(in)>=0))
      return false;
    value = static_cast<std::size_t>(PyLong_AS_LONG(in));
    return true;
%#else
    if (!(PyInt_Check(in) && PyInt_AS_LONG(in)>=0))
      return false;
    value = static_cast<std::size_t>(PyInt_AS_LONG(in));
    return true;
%#endif
  }
}

//-----------------------------------------------------------------------------
// Out typemap (std::size_t)
//-----------------------------------------------------------------------------
%typemap(out) std::size_t
{
  // Typemap std::size_t
%#if PY_MAJOR_VERSION >= 3
    $result = PyLong_FromLong(static_cast< long >($1));
%#else
    $result = PyInt_FromLong(static_cast< long >($1));
%#endif
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
del UFC_VERSION, UFC_VERSION_MAJOR, UFC_VERSION_MINOR, UFC_VERSION_MAINTENANCE, UFC_VERSION_RELEASE

"""Code for adding swig version to ufc extension module."""
tmp = hex(ufc_swigversion())
__swigversion__ = "%d.%d.%d"%(tuple(map(int, [tmp[-5], tmp[-3], tmp[-2:]])))

del tmp, ufc_swigversion
%}
