//
// This code complies with UFC version 1.0, and is generated with SyFi version 0.4.0.
//
// http://www.fenics.org/syfi/
// http://www.fenics.org/ufc/
//


#ifndef __form__stiffness_form__Lagrange_1_2D_H
#define __form__stiffness_form__Lagrange_1_2D_H

#include <stdexcept>
#include <math.h>
#include <ufc.h>
#include <pycc/Functions/Ptv.h>
#include <pycc/Functions/Ptv_tools.h>
#include <pycc/Functions/Dof_Ptv.h>
#include <pycc/Functions/OrderedPtvSet.h>
#include <pycc/Functions/Dof_OrderedPtvSet.h>
#include "dof_map_Lagrange_1_2D.h"
#include "fe_Lagrange_1_2D.h"


namespace pycc
{

/// This class defines the interface for the assembly of the global
/// tensor corresponding to a form with r + n arguments, that is, a
/// mapping
///
///     a : V1 x V2 x ... Vr x W1 x W2 x ... x Wn -> R
///
/// with arguments v1, v2, ..., vr, w1, w2, ..., wn. The rank r
/// global tensor A is defined by
///
///     A = a(V1, V2, ..., Vr, w1, w2, ..., wn),
///
/// where each argument Vj represents the application to the
/// sequence of basis functions of Vj and w1, w2, ..., wn are given
/// fixed functions (coefficients).

class form__stiffness_form__Lagrange_1_2D: public ufc::form
{  
public:

  /// Constructor
  form__stiffness_form__Lagrange_1_2D();

  /// Destructor
  virtual ~form__stiffness_form__Lagrange_1_2D();

  /// Return a string identifying the form
  virtual const char* signature() const;

  /// Return the rank of the global tensor (r)
  virtual unsigned int rank() const;

  /// Return the number of coefficients (n)
  virtual unsigned int num_coefficients() const;

  /// Return the number of cell integrals
  virtual unsigned int num_cell_integrals() const;
  
  /// Return the number of exterior facet integrals
  virtual unsigned int num_exterior_facet_integrals() const;
  
  /// Return the number of interior facet integrals
  virtual unsigned int num_interior_facet_integrals() const;

  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(unsigned int i) const;

  /// Create a new dof map for argument function i
  virtual ufc::dof_map* create_dof_map(unsigned int i) const;

  /// Create a new cell integral on sub domain i
  virtual ufc::cell_integral* create_cell_integral(unsigned int i) const;

  /// Create a new exterior facet integral on sub domain i
  virtual ufc::exterior_facet_integral* 
    create_exterior_facet_integral(unsigned int i) const;

  /// Create a new interior facet integral on sub domain i
  virtual ufc::interior_facet_integral* 
    create_interior_facet_integral(unsigned int i) const;

};


} // namespace

#endif
