// This is UFC (Unified Form-assembly Code) specification v. 0.0.
// This code has been released into the public domain.
//
// Simula Research Laboratory, 2006.

// FIXME: --> The FEniCS project?
// FIXME: --> Finite Element Center?

#ifndef __UFC_H
#define __UFC_H

namespace ufc
{

  /// This class defines the data structure for a finite element mesh. 

  class mesh
  {
  public:
    mesh(): num_entities(0) {}

    /// Space dimension and length of num_entities
    unsigned int space_dimension;
    
    /// Array of the number of entities for each topological dimension
    unsigned int* num_entities;

  };
  
  /// This class defines the data structure for a cell in a
  /// finite element mesh.

  class cell
  {
  public:
    cell(): entities(0), coordinates(0) {}

    /// Array of global indices for mesh entities contained in the cell
    unsigned int** entities;

    /// Array of coordinates for the vertices of the cell
    double** coordinates;
    
  };

  /// This class defines the interface for the local-to-global mapping
  /// of nodes (degrees of freedom) for a finite element space.

  class node_map
  {
  public:
    /// Destructor
    virtual ~node_map() {}

    /// Return dimension of the local finite element space
    virtual unsigned int space_dimension() const = 0;

    /// Tabulate local-to-global mapping of nodes (degrees of freedom)
    virtual void tabulate(unsigned int *nodes, const mesh& m, const cell& c) const = 0;
 
  };


  /// This class defines a generic function f: Input -> Output
  // TODO: Define Input and Output statically? Which type?

  template<class Input, class Output>
  class function
  {
  public:
    /// Destructor
    virtual ~function() {}

    /// Evaluating
    virtual Output eval(Input) = 0;
  };


  /// This class defines a mapping of a generic function<Input, Output>
  /// onto a local finite element space.

  template<class Input, class Output>
  class projection
  {
  public:
    /// Destructor
    virtual ~projection() {}

    /// For a particular finite element, computes the nodal functional of f:
    /// dofs[i] = nu_i(f)  for each local dof i
    virtual void project(function<Input,Output> & f, const cell & c, double *dofs) = 0;
  };


  /// This class defines the interface for the computation of the
  /// element tensor corresponding to a form with r + n arguments,
  /// that is, a mapping
  ///
  ///     a : V1 x V2 x ... Vr x W1 x W2 x ... x Wn -> R
  ///
  /// with arguments v1, v2, ..., vr, w1, w2, ..., wn. The rank r
  /// element tensor AK is the local contribution from an element K
  /// to the rank r global tensor A defined by
  ///
  ///     Ai = a(V1, V2, ..., Vr, w1, w2, ..., wn),
  ///
  /// where each argument Vj represents the application to the sequence
  /// of nodal basis functions of Vj and w1, w2, ..., wn are given fixed
  /// functions (coefficients).
  
  class element_tensor
  {
  public:

    /// Constructor
    element_tensor() : node_maps(0) {}

    /// Destructor
    virtual ~element_tensor() {}

    /// Return the rank of the tensor (r)
    virtual unsigned int rank() const = 0;

    /// Return the number of coefficients (n)
    virtual unsigned int num_coefficients() const = 0;

    /// Tabulate element tensor
    virtual void tabulate(double* A, const double** w, const cell& c) const = 0;
 
    /// Array of r + n node maps for argument functions and coefficients
    node_map** node_maps;
  
  };

}

#endif
