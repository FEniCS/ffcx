// This is UFC (Unified Form-assembly Code) v. 0.0.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenics.org/) 2006.

#ifndef __UFC_H
#define __UFC_H

// FIXME: Should we have constructors and destructors for each class?

namespace ufc
{

  /// This class defines the data structure for a finite element mesh. 

  class mesh
  {
  public:

    /// Constructor
    mesh(): num_entities(0) {}

    // FIXME: Do we need this? Length of num_entities is the topological
    // FIXME: dimension + 1.
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

    /// Constructor
    cell(): entities(0), coordinates(0) {}

    /// Array of global indices for mesh entities contained in the cell
    unsigned int** entities;

    /// Array of coordinates for the vertices of the cell
    double** coordinates;
    
  };

  /// This class defines the interface for a tensor-valued function.

  class function
  {
  public:

    /// Destructor
    virtual ~function() {}

    /// Evaluate component i = (i[0], i[1], ...) at the point x = (x[0], x[1], ...)
    virtual double eval(const int* i, const double* x) const;

  };

  /// This class defines the interface for a finite element.

  class finite_element
  {
  public:

    /// Destructor
    virtual ~finite_element() {}

    /// Return dimension of the finite element function space
    virtual unsigned int space_dimension() const = 0;

    /// Evaluate basis function n at given point in cell
    virtual double evaluate_basis(unsigned int n, double* x, const cell& c) const = 0;
    
    /// Evaluate node n at given function
    virtual double evaluate_node(unsigned int n, const function& f, const cell& c) const = 0;

    /// Compute local-to-global mapping of nodes (degrees of freedom)
    virtual void tabulate_nodes(unsigned int *nodes, const mesh& m, const cell& c) const = 0;

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
    element_tensor() : finite_elements(0) {}

    /// Destructor
    virtual ~element_tensor() {}

    /// Return the rank of the element tensor (r)
    virtual unsigned int rank() const = 0;

    /// Return the number of coefficients (n)
    virtual unsigned int num_coefficients() const = 0;

    /// Return true iff there is a contribution from the interior
    virtual bool interior_contribution() const = 0;

    /// Return true iff there is a contribution from the boundary
    virtual bool boundary_contribution() const = 0;

    /// Tabulate interior contribution to element tensor
    virtual void tabulate_interior(double* A, const double** w, const cell& c) const = 0;

    /// Tabulate boundary contribution to element tensor
    virtual void tabulate_boundary(double* A, const double** w, const cell& c, unsigned int facet) const = 0;
 
    /// Array of r + n finite elements for argument function spaces and coefficients
    finite_element** finite_elements;
  
  };

}

#endif
