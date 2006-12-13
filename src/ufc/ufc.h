// This is UFC (Unified Form-assembly Code) v. 1.0-rc1
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
    
  /// This class defines the interface for a local-to-global mapping of
  /// degrees of freedom (dofs).

  class dof_map
  {
  public:

    /// Constructor
    dof_map() {}

    /// Destructor
    virtual ~dof_map() {}
    
    /// Return a string identifying the dof map
    virtual const char* signature() const = 0;
    
    /// Return true iff mesh entities of topological dimension d are needed
    virtual bool needs_mesh_entities(unsigned int d) const = 0;

    /// Initialize dof map for a new mesh and return true iff
    /// the finite element should be initialized for each cell
    virtual bool init_mesh(const mesh& mesh) = 0;

    /// Initialize dof map for given cell
    virtual void init_cell(const mesh& mesh, const cell& cell) = 0;

    /// Return the dimension of the global finite element function space.
    /// This is only valid after initialization.
    virtual unsigned int global_dimension() const = 0;
    
    /// Return the dimension of the local finite element function space
    virtual unsigned int local_dimension() const = 0;
    
    /// Tabulate the local-to-global mapping of dofs
    virtual void tabulate_dofs(unsigned int *dofs, const mesh& m, const cell& c) const = 0;
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
    virtual const char* signature() const = 0;

    /// Return the dimension of the finite element function space
    virtual unsigned int space_dimension() const = 0;

    /// Return the rank of the value space
    virtual unsigned int value_rank() const = 0;

    /// Return the dimension of the value space for axis i
    virtual unsigned int value_dimension(unsigned int i) const = 0;

    /// Evaluate basis function i at the point x = (x[0], x[1], ...) in the cell
    virtual void evaluate_basis(double* values, const double* x, unsigned int i, const cell& c) const = 0;
    
    /// Evaluate linear functional for dof i on the function f
    virtual double evaluate_dof(unsigned int i, const function& f, const cell& c) const = 0;
    
    /// Interpolate vertex values from dof values
    virtual void interpolate_vertex_values(double* vertex_values, const double* dof_values) const = 0;
  
    /// Return the number of sub elements (for a mixed finite element)
    virtual unsigned int num_sub_elements(unsigned int i) const = 0;

    /// Return sub element i (for a mixed finite element)
    virtual const finite_element& sub_element(unsigned int i) const = 0;
    
  };

  /// This class defines the interface for the tabulation of the cell
  /// tensor corresponding to the local contribution to a form from
  /// the integral over a cell.

  class cell_integral
  {
  public:

    /// Constructor
    cell_integral() {}

    /// Destructor
    virtual ~cell_integral() {}

    /// Tabulate the tensor for the contribution from a local cell
    virtual void tabulate_tensor(double* A, const double * const * w, const cell& c) const = 0;

  };

  /// This class defines the interface for the tabulation of the
  /// exterior facet tensor corresponding to the local contribution
  /// to a form from the integral over an exterior facet.

  class exterior_facet_integral
  {
  public:

    /// Constructor
    exterior_facet_integral() {}

    /// Destructor
    virtual ~exterior_facet_integral() {}

    /// Tabulate the tensor for the contribution from a local exterior facet
    virtual void tabulate_tensor(double* A, const double * const * w, const cell& c, unsigned int facet) const = 0;

  };

  /// This class defines the interface for the tabulation of the
  /// interior facet tensor corresponding to the local contribution
  /// to a form from the integral over an interior facet.

  class interior_facet_integral
  {
  public:

    /// Constructor
    interior_facet_integral() {}

    /// Destructor
    virtual ~interior_facet_integral() {}

    /// Tabulate the tensor for the contribution from a local interior facet
    virtual void tabulate_tensor(double* A, const double * const * w, const cell& c0, const cell& c1, unsigned int facet0, unsigned int facet1) const = 0;

  };

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
  /// where each argument Vj represents the application to the sequence
  /// of basis functions of Vj and w1, w2, ..., wn are given fixed
  /// functions (coefficients).
  
  class form
  {
  public:

    /// Constructor
    form() {}

    /// Destructor
    virtual ~form() {}

    /// Return a string identifying the form
    virtual const char* signature() const = 0;

    /// Return the rank of the element tensor (r)
    virtual unsigned int rank() const = 0;

    /// Return the number of coefficients (n)
    virtual unsigned int num_coefficients() const = 0;

    /// Create cell integral, returning 0 if the contribution is zero (caller responsible for deletion)
    virtual cell_integral* create_cell_integral() const = 0;
    
    /// Create interior facet integral, returning 0 if the contribution is zero (caller responsible for deletion)
    virtual interior_facet_integral* create_interior_facet_integral() const = 0;

    /// Create exterior facet integral, returning 0 if the contribution is zero (caller responsible for deletion)
    virtual exterior_facet_integral* create_exterior_facet_integral() const = 0;
    
    /// Create a dof map for argument function i (caller responsible for deletion)
    virtual dof_map* create_dof_map(unsigned int i) const = 0;

    /// Create a finite element for argument function i (caller responsible for deletion)
    virtual finite_element* create_finite_element(unsigned int i) const = 0;

  };

}

#endif
