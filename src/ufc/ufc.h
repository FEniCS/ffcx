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
    mesh(): topological_dimension(0), geometric_dimension(0), num_entities(0) {}

    /// Destructor
    virtual ~mesh() {}

    /// Topological dimension of the mesh
    unsigned int topological_dimension;
    
    /// Geometric dimension of the mesh
    unsigned int geometric_dimension;

    /// Array of the global number of entities of each topological dimension
    unsigned int* num_entities;

  };
  
  /// This class defines the data structure for a cell in a mesh.
  
  class cell
  {
  public:

    /// Constructor
    cell(): entity_indices(0), coordinates(0) {}

    /// Destructor
    virtual ~cell() {}

    /// Array of global indices for the mesh entities of the cell
    unsigned int** entity_indices;

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

    /// Evaluate the function at the point x = (x[0], x[1], ...) in the cell
    virtual void evaluate(double* values, const double* x, const cell& c) const = 0;

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

    /// Return true iff indices of mesh entities with topological dimension d are needed
    virtual bool needs_mesh_entities(unsigned int d) const = 0;

    /// Initialize finite element for a new mesh and return true iff
    /// the finite element should be initialized on each cell
    virtual bool init(const mesh& mesh) = 0;

    /// Initialize finite element on given cell
    virtual void init(const mesh& mesh, const cell& cell) = 0;

    /// Return the dimension of the global finite element function space.
    /// This is only valid after calling init(mesh) for the relevant mesh.
    virtual unsigned int global_dimension() const = 0;

    /// Return the dimension of the local finite element function space
    virtual unsigned int local_dimension() const = 0;

    /// Return the rank of the value space
    virtual unsigned int value_rank() const = 0;

    /// Return the dimension of value space for axis i
    virtual unsigned int value_dimension(unsigned int i) const = 0;

    /// Return the number of sub elements (for a mixed finite element)
    virtual unsigned int num_sub_elements(unsigned int i) const = 0;

    /// Return sub element i (for a mixed finite element)
    virtual const finite_element& sub_element(unsigned int i) const = 0;
    
    /// Evaluate basis function i at the point x = (x[0], x[1], ...) in the cell
    virtual void evaluate_basis(double* values, const double* x, unsigned int i, const cell& c) const = 0;
    
    /// Evaluate node i on the function f
    virtual double evaluate_node(unsigned int n, const function& f, const cell& c) const = 0;
    
    /// Interpolate vertex values from nodal values
    virtual void interpolate_vertex_values(double* vertex_values, const double* nodal_values) const = 0;

    /// Tabulate the local-to-global mapping of nodes (degrees of freedom)
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

    /// Return a string identifying the form
    virtual const char* description() const = 0;

    /// Return the rank of the element tensor (r)
    virtual unsigned int rank() const = 0;

    /// Return the number of coefficients (n)
    virtual unsigned int num_coefficients() const = 0;

    /// Return true iff form has a contribution from the interior
    virtual bool has_interior_contribution() const = 0;

    /// Return true iff form has a contribution from the boundary
    virtual bool has_boundary_contribution() const = 0;

    /// Tabulate the interior contribution to the element tensor
    virtual void tabulate_interior(double* A, const double * const * w, const cell& c) const = 0;
    
    /// Tabulate the boundary contribution to the element tensor
    virtual void tabulate_boundary(double* A, const double * const * w, const cell& c, unsigned int facet) const = 0;
 
    /// Array of r + n finite elements for argument function spaces and coefficients
    finite_element** finite_elements;
  
  };

}

#endif
