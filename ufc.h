// This is UFC (Unified Form-assembly Code) v. 0.0.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenics.org/) 2006.

#ifndef __UFC_H
#define __UFC_H

namespace ufc
{

  /// This class defines the data structure for a finite element mesh.

  class mesh
  {
  public:

    /// Constructor
    mesh(): num_entities(0) {}

    /// Destructor
    ~mesh() {}

    /// Array of the number of entities of each topological dimension
    unsigned int* num_entities;

  };
  
  /// This class defines the data structure for a cell in a finite
  /// element mesh.

  class cell
  {
  public:

    /// Constructor
    cell(): entities(0), coordinates(0) {}

    /// Destructor
    ~cell() {}

    /// Array of global indices for the mesh entities of the cell
    unsigned int** entities;

    /// Array of coordinates for the vertices of the cell
    double** coordinates;
    
  };

  /// This class defines the interface for a general tensor-valued function.

  class function
  {
  public:

    /// Constructor
    function() {}

    /// Destructor
    virtual ~function() {}

    /// Evaluate component i = (i[0], i[1], ...) at the point x = (x[0], x[1], ...)
    virtual double eval(const unsigned int* i, const double* x) const = 0;

  };

  /// This class defines the interface for a finite element.

  class finite_element
  {
  public:

    /// Constructor
    finite_element() {}

    /// Destructor
    virtual ~finite_element() {}

    /// Return a string identifying the finite element
    virtual const char* description() const = 0;

    /// Return dimension of the finite element function space
    virtual unsigned int space_dimension() const = 0;

    /// Return rank of value space
    virtual unsigned int value_rank() const = 0;

    /// Return dimension of value space for axis i
    virtual unsigned int value_dimension(unsigned int i) const = 0;

    /// Return number of sub elements (for a mixed finite element)
    virtual unsigned int num_sub_elements(unsigned int i) const = 0;

    /// Return sub element i (for a mixed finite element)
    virtual const finite_element& sub_element(unsigned int i) const = 0;
    
    /// Evaluate basis function n at a given point in the cell
    virtual double evaluate_basis(unsigned int n, double* x, const cell& c) const = 0;
    
    /// Evaluate node n at given function
    virtual double evaluate_node(unsigned int n, const function& f, const cell& c) const = 0;

    /// Tabulate local-to-global mapping of nodes (degrees of freedom)
    virtual void tabulate_nodes(unsigned int *nodes, const mesh& m, const cell& c) const = 0;

    /// Tabulate values at the vertices of the cell
    virtual void tabulate_vertex_values(double* vertex_values, const double* nodal_values) const = 0;
    
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

    /// Return a string identifying the form
    virtual const char* description() const = 0;

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
