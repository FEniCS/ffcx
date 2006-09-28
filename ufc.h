// Copyright (C) 2006 Anders Logg, Kent-Andre Mardal, Ola Skavhaug.
// Licensed under the GNU GPL Version 2.

#ifndef __UFC_H
#define __UFC_H

namespace ufc
{

  /// This class defines the interface for a finite element generating
  /// the function space for an  argument function or coefficient in a
  /// multilinear form.   It does not define a  complete interface for
  /// the standard  concept of  a finite element,  but only  a minimal
  /// interface  for the  functionality that  may be  requested during
  /// assembly of the discrete representation of the multilinear form.

  class finite_element
  {
  public:
    
    /// Return dimension of the local finite element space
    virtual unsigned int space_dimension() const = 0;

    /// Return dimension of the local domain
    virtual unsigned int shape_dimension() const = 0;

    /// Compute local-to-global mapping of nodes (degrees of freedom)
    virtual void compute_node_map() const = 0;
    
  };

  /// This class defines the interface for a multilinear form, that is
  /// a mapping
  ///
  ///     a : V1 x V2 x ... Vr x W1 x W2 x ... x Wn -> R
  ///
  ///     a = a(v1, v2, ..., vr; w1, w2, ..., wn)
  ///
  /// with arguments v1, v2, ..., vr and coefficients w1, w2, ..., wn.
  ///
  /// The interface is limited to the minimal set of functions that
  /// are needed to assemble the discrete representation of the form.
  
  class multilinear_form
  {
  public:

    /// Constructor
    multilinear_form() : finite_elements(0) {}

    /// Destructor
    virtual ~multilinear_form() {}
    
    /// Return number of argument functions (rank, arity)
    virtual unsigned int num_arguments() const = 0;

    /// Return number of coefficients
    virtual unsigned int num_coefficients() const = 0;

    /// Return whether there is a contribution from the interior
    virtual bool interior_contribution() const = 0;

    /// Return whether there is a contribution from the boundary
    virtual bool boundary_contribution() const = 0;

    /// Compute interior contribution to element tensor
    virtual void compute_element_tensor_interior(double* A) const = 0;
    
    /// Compute boundary contribution to element tensor
    virtual void compute_element_tensor_boundary(double* A) const = 0;

    /// Array of finite elements for argument functions and coefficients
    finite_element** finite_elements;

  };

}

#endif
