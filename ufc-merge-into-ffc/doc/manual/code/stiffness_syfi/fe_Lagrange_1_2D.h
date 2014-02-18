//
// This code complies with UFC version 1.0, and is generated with SyFi version 0.4.0.
//
// http://www.fenics.org/syfi/
// http://www.fenics.org/ufc/
//


#ifndef __fe_Lagrange_1_2D_H
#define __fe_Lagrange_1_2D_H

#include <stdexcept>
#include <math.h>
#include <ufc.h>
#include <pycc/Functions/Ptv.h>
#include <pycc/Functions/Ptv_tools.h>
#include <pycc/Functions/Dof_Ptv.h>
#include <pycc/Functions/OrderedPtvSet.h>
#include <pycc/Functions/Dof_OrderedPtvSet.h>



namespace pycc
{

/// This class defines the interface for a finite element.

class fe_Lagrange_1_2D: public ufc::finite_element
{  
public:

  /// Constructor
  fe_Lagrange_1_2D();

  /// Destructor
  virtual ~fe_Lagrange_1_2D();

  /// Return a string identifying the finite element
  virtual const char* signature() const;

  /// Return the cell shape
  virtual ufc::shape cell_shape() const;

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const;

  /// Return the rank of the value space
  virtual unsigned int value_rank() const;

  /// Return the dimension of the value space for axis i
  virtual unsigned int value_dimension(unsigned int i) const;

  /// Evaluate basis function i at given point in cell
  virtual void evaluate_basis(unsigned int i,
                              double* values,
                              const double* coordinates,
                              const ufc::cell& c) const;

  /// Evaluate order n derivatives of basis function i at given point in cell
  virtual void evaluate_basis_derivatives(unsigned int i,
                                          unsigned int n,
                                          double* values,
                                          const double* coordinates,
                                          const ufc::cell& c) const;

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const;

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const;

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const;

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const;

};


} // namespace

#endif
