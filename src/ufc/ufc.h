// This is UFC (Unified Form-assembly Code) v. 1.4.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2006-2010.

#ifndef __UFC_H
#define __UFC_H

#define UFC_VERSION_MAJOR 1
#define UFC_VERSION_MINOR 4
#define UFC_VERSION_MAINTENANCE 2

#include <stdexcept>

const char UFC_VERSION[] = "1.4";

namespace ufc
{

  /// Valid cell shapes
  enum shape {interval, triangle, quadrilateral, tetrahedron, hexahedron};

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
    cell(): cell_shape(interval),
            topological_dimension(0), geometric_dimension(0),
            entity_indices(0), coordinates(0), index(0), local_facet(-1),
            mesh_identifier(-1) {}

    /// Destructor
    virtual ~cell() {}

    /// Shape of the cell
    shape cell_shape;

    /// Topological dimension of the mesh
    unsigned int topological_dimension;

    /// Geometric dimension of the mesh
    unsigned int geometric_dimension;

    /// Array of global indices for the mesh entities of the cell
    unsigned int** entity_indices;

    /// Array of coordinates for the vertices of the cell
    double** coordinates;

    /// Cell index (short-cut for entity_indices[topological_dimension][0])
    unsigned int index;

    /// Local facet index
    int local_facet;

    /// Unique mesh identifier
    int mesh_identifier;

  };

  /// This class defines the interface for a general tensor-valued function.

  class function
  {
  public:

    /// Constructor
    function() {}

    /// Destructor
    virtual ~function() {}

    /// Evaluate function at given point in cell
    virtual void evaluate(double* values,
                          const double* coordinates,
                          const cell& c) const = 0;

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

    /// Return the cell shape
    virtual shape cell_shape() const = 0;

    /// Return the dimension of the finite element function space
    virtual unsigned int space_dimension() const = 0;

    /// Return the rank of the value space
    virtual unsigned int value_rank() const = 0;

    /// Return the dimension of the value space for axis i
    virtual unsigned int value_dimension(unsigned int i) const = 0;

    /// Evaluate basis function i at given point in cell
    virtual void evaluate_basis(unsigned int i,
                                double* values,
                                const double* coordinates,
                                const cell& c) const = 0;

    /// Evaluate all basis functions at given point in cell
    virtual void evaluate_basis_all(double* values,
                                    const double* coordinates,
                                    const cell& c) const = 0;

    /// Evaluate order n derivatives of basis function i at given point in cell
    virtual void evaluate_basis_derivatives(unsigned int i,
                                            unsigned int n,
                                            double* values,
                                            const double* coordinates,
                                            const cell& c) const = 0;

    /// Evaluate order n derivatives of all basis functions at given point in cell
    virtual void evaluate_basis_derivatives_all(unsigned int n,
                                                double* values,
                                                const double* coordinates,
                                                const cell& c) const = 0;

    /// Evaluate linear functional for dof i on the function f
    virtual double evaluate_dof(unsigned int i,
                                const function& f,
                                const cell& c) const = 0;

    /// Evaluate linear functionals for all dofs on the function f
    virtual void evaluate_dofs(double* values,
                               const function& f,
                               const cell& c) const = 0;

    /// Interpolate vertex values from dof values
    virtual void interpolate_vertex_values(double* vertex_values,
                                           const double* dof_values,
                                           const cell& c) const = 0;

    /// Return the number of sub elements (for a mixed element)
    virtual unsigned int num_sub_elements() const = 0;

    /// Create a new finite element for sub element i (for a mixed element)
    virtual finite_element* create_sub_element(unsigned int i) const = 0;

#ifdef UFC_DEV

    /// Create a new class instance
    virtual finite_element* create() const
    {
      throw std::runtime_error("finite_element::create not implemented (introduced in UFC v1.6).");
    }

#endif

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

    /// Initialize dof map for mesh (return true iff init_cell() is needed)
    virtual bool init_mesh(const mesh& mesh) = 0;

    /// Initialize dof map for given cell
    virtual void init_cell(const mesh& m,
                           const cell& c) = 0;

    /// Finish initialization of dof map for cells
    virtual void init_cell_finalize() = 0;

    /// Return the dimension of the global finite element function space
    virtual unsigned int global_dimension() const = 0;

    /// Return the dimension of the local finite element function space for a cell
    virtual unsigned int local_dimension(const cell& c) const = 0;

    /// Return the maximum dimension of the local finite element function space
    virtual unsigned int max_local_dimension() const = 0;

    // Return the geometric dimension of the coordinates this dof map provides
    virtual unsigned int geometric_dimension() const = 0;

    /// Return the number of dofs on each cell facet
    virtual unsigned int num_facet_dofs() const = 0;

    /// Return the number of dofs associated with each cell entity of dimension d
    virtual unsigned int num_entity_dofs(unsigned int d) const = 0;

    /// Tabulate the local-to-global mapping of dofs on a cell
    virtual void tabulate_dofs(unsigned int* dofs,
                               const mesh& m,
                               const cell& c) const = 0;

    /// Tabulate the local-to-local mapping from facet dofs to cell dofs
    virtual void tabulate_facet_dofs(unsigned int* dofs,
                                     unsigned int facet) const = 0;

    /// Tabulate the local-to-local mapping of dofs on entity (d, i)
    virtual void tabulate_entity_dofs(unsigned int* dofs,
                                      unsigned int d, unsigned int i) const = 0;

    /// Tabulate the coordinates of all dofs on a cell
    virtual void tabulate_coordinates(double** coordinates,
                                      const cell& c) const = 0;

    /// Return the number of sub dof maps (for a mixed element)
    virtual unsigned int num_sub_dof_maps() const = 0;

    /// Create a new dof_map for sub dof map i (for a mixed element)
    virtual dof_map* create_sub_dof_map(unsigned int i) const = 0;

#ifdef UFC_DEV

    /// Create a new class instance
    virtual dof_map* create() const
    {
      throw std::runtime_error("dof_map::create not implemented (introduced in UFC v1.6).");
    }

#endif

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
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const cell& c) const = 0;

  };

  /// This class defines the interface for the tabulation of the
  /// exterior facet tensor corresponding to the local contribution to
  /// a form from the integral over an exterior facet.

  class exterior_facet_integral
  {
  public:

    /// Constructor
    exterior_facet_integral() {}

    /// Destructor
    virtual ~exterior_facet_integral() {}

    /// Tabulate the tensor for the contribution from a local exterior facet
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const cell& c,
                                 unsigned int facet) const = 0;

  };

  /// This class defines the interface for the tabulation of the
  /// interior facet tensor corresponding to the local contribution to
  /// a form from the integral over an interior facet.

  class interior_facet_integral
  {
  public:

    /// Constructor
    interior_facet_integral() {}

    /// Destructor
    virtual ~interior_facet_integral() {}

    /// Tabulate the tensor for the contribution from a local interior facet
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const cell& c0,
                                 const cell& c1,
                                 unsigned int facet0,
                                 unsigned int facet1) const = 0;

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
  /// where each argument Vj represents the application to the
  /// sequence of basis functions of Vj and w1, w2, ..., wn are given
  /// fixed functions (coefficients).

  class form
  {
  public:

    /// Constructor
    form() {}

    /// Destructor
    virtual ~form() {}

    /// Return a string identifying the form
    virtual const char* signature() const = 0;

    /// Return the rank of the global tensor (r)
    virtual unsigned int rank() const = 0;

    /// Return the number of coefficients (n)
    virtual unsigned int num_coefficients() const = 0;

    /// Return the number of cell integrals
    virtual unsigned int num_cell_integrals() const = 0;

    /// Return the number of exterior facet integrals
    virtual unsigned int num_exterior_facet_integrals() const = 0;

    /// Return the number of interior facet integrals
    virtual unsigned int num_interior_facet_integrals() const = 0;

    /// Create a new finite element for argument function i
    virtual finite_element* create_finite_element(unsigned int i) const = 0;

    /// Create a new dof map for argument function i
    virtual dof_map* create_dof_map(unsigned int i) const = 0;

    /// Create a new cell integral on sub domain i
    virtual cell_integral* create_cell_integral(unsigned int i) const = 0;

    /// Create a new exterior facet integral on sub domain i
    virtual exterior_facet_integral*
    create_exterior_facet_integral(unsigned int i) const = 0;

    /// Create a new interior facet integral on sub domain i
    virtual interior_facet_integral*
    create_interior_facet_integral(unsigned int i) const = 0;

  };

}

#endif
