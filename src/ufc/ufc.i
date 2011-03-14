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
%ignore ufc::cell::coordinates;

%include "ufc.h"

// Declare which classes should be stored using shared_ptr
%shared_ptr(ufc::cell_integral)
%shared_ptr(ufc::dofmap)
%shared_ptr(ufc::finite_element)
%shared_ptr(ufc::function)
%shared_ptr(ufc::form)
%shared_ptr(ufc::exterior_facet_integral)
%shared_ptr(ufc::interior_facet_integral)

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
