%module ufc

%{
#include <ufc.h>
// If using std::tr1::shared_ptr comment out this line
#include <boost/shared_ptr.hpp>
// Un comment these lines to use std::tr1, only works with swig version  >= 1.3.37
//#include <tr1/memory>
%}


// Handle shared_ptr only available for swig version >= 1.3.34
#if SWIG_VERSION >= 0x010334
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

%include ufc.h

#endif
