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
%include "boost_shared_ptr.i"

SWIG_SHARED_PTR(form,ufc::form)
SWIG_SHARED_PTR(finite_element,ufc::finite_element)
SWIG_SHARED_PTR(dof_map,ufc::dof_map)
SWIG_SHARED_PTR(cell_integral,ufc::cell_integral)
SWIG_SHARED_PTR(exterior_facet_integral,ufc::exterior_facet_integral)
SWIG_SHARED_PTR(interior_facet_integral,ufc::interior_facet_integral)

%include "ufc.h"

// Include code to generate a __swigversion__ attribute to the cpp module
// Add prefix to avoid naming problems with other modules

%inline %{
int ufc_swigversion() { return SWIGVERSION; }
%}

int ufc_swigversion();

%pythoncode %{
"""Preliminary code for adding swig version to cpp module. Someone (tm) finish
this.
"""
tmp = hex(ufc_swigversion())
__swigversion__ = ".".join([tmp[-5],tmp[-3],tmp[-2:]])

del tmp, ufc_swigversion
%}
