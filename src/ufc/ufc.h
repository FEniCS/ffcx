// This is UFC (Unified Form-assembly Code) v. 2.3.0.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2006-2014.

#ifndef __UFC_H
#define __UFC_H

#define UFC_VERSION_MAJOR 2
#define UFC_VERSION_MINOR 3
#define UFC_VERSION_MAINTENANCE 0

#include <vector>
#include <cstddef>
#include <stdexcept>

#include <ufc_geometry.h>

const char UFC_VERSION[] = "2.3.0";

namespace ufc
{

  /// Valid cell shapes
  enum shape {interval, triangle, quadrilateral, tetrahedron, hexahedron};

  /// This class defines the interface for cell topology data.
  class cell_topology
  {
  public:

    /// Destructor
    virtual ~cell_topology() {}

    /// Return array of global entity indices for topological dimension d
    //virtual const std::size_t* entity_indices(std::size_t d) const;

  };

  /// This class defines the interface for cell geometry data.
  class cell_geometry
  {
  public:

    /// Destructor
    virtual ~cell_geometry() {}

    /// Get vertex coordinates
    //virtual void get_vertex_coordinates(const double* x[]) const;

    //virtual std::vector<std::vector<double> >& void vertex_coordinates() const;

  };

  /// This class defines the data structure for a cell in a mesh.

  class cell
  {
  public:

    /// Constructor
    cell(): cell_shape(interval),
            topological_dimension(0), geometric_dimension(0),
            index(0), local_facet(-1), mesh_identifier(-1) {}

    /// Destructor
    virtual ~cell() {}

    /// Shape of the cell
    shape cell_shape;

    /// Topological dimension of the mesh
    std::size_t topological_dimension;

    /// Geometric dimension of the mesh
    std::size_t geometric_dimension;

    /// Array of global indices for the mesh entities of the cell
    std::vector<std::vector<std::size_t> > entity_indices;

    /// Cell index (short-cut for entity_indices[topological_dimension][0])
    std::size_t index;

    /// Local facet index
    int local_facet;

    /// Cell orientation
    int orientation;

    /// Unique mesh identifier
    int mesh_identifier;

  };

  /// This class defines the interface for a general tensor-valued function.

  class function
  {
  public:

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

    /// Destructor
    virtual ~finite_element() {}

    /// Return a string identifying the finite element
    virtual const char* signature() const = 0;

    /// Return the cell shape
    virtual shape cell_shape() const = 0;

    /// Return the topological dimension of the cell shape
    virtual std::size_t topological_dimension() const = 0;

    /// Return the geometric dimension of the cell shape
    virtual std::size_t geometric_dimension() const = 0;

    /// Return the dimension of the finite element function space
    virtual std::size_t space_dimension() const = 0;

    /// Return the rank of the value space
    virtual std::size_t value_rank() const = 0;

    /// Return the dimension of the value space for axis i
    virtual std::size_t value_dimension(std::size_t i) const = 0;

    /// Evaluate basis function i at given point x in cell
    virtual void evaluate_basis(std::size_t i,
                                double* values,
                                const double* x,
                                const double* vertex_coordinates,
                                int cell_orientation) const = 0;

    /// Evaluate all basis functions at given point x in cell
    virtual void evaluate_basis_all(double* values,
                                    const double* x,
                                    const double* vertex_coordinates,
                                    int cell_orientation) const = 0;

    /// Evaluate order n derivatives of basis function i at given point x in cell
    virtual void evaluate_basis_derivatives(std::size_t i,
                                            std::size_t n,
                                            double* values,
                                            const double* x,
                                            const double* vertex_coordinates,
                                            int cell_orientation) const = 0;

    /// Evaluate order n derivatives of all basis functions at given point x in cell
    virtual void evaluate_basis_derivatives_all(std::size_t n,
                                                double* values,
                                                const double* x,
                                                const double* vertex_coordinates,
                                                int cell_orientation) const = 0;

    // FIXME: cell argument only included here so we can pass it to the eval function...

    /// Evaluate linear functional for dof i on the function f
    virtual double evaluate_dof(std::size_t i,
                                const function& f,
                                const double* vertex_coordinates,
                                int cell_orientation,
                                const cell& c) const = 0;

    /// Evaluate linear functionals for all dofs on the function f
    virtual void evaluate_dofs(double* values,
                               const function& f,
                               const double* vertex_coordinates,
                               int cell_orientation,
                               const cell& c) const = 0;

    /// Interpolate vertex values from dof values
    virtual void interpolate_vertex_values(double* vertex_values,
                                           const double* dof_values,
                                           const double* vertex_coordinates,
                                           int cell_orientation,
                                           const cell& c) const = 0;

    /// Map coordinate xhat from reference cell to coordinate x in cell
    virtual void map_from_reference_cell(double* x,
                                         const double* xhat,
                                         const cell& c) const = 0;

    /// Map from coordinate x in cell to coordinate xhat in reference cell
    virtual void map_to_reference_cell(double* xhat,
                                       const double* x,
                                       const cell& c) const = 0;

    /// Return the number of sub elements (for a mixed element)
    virtual std::size_t num_sub_elements() const = 0;

    /// Create a new finite element for sub element i (for a mixed element)
    virtual finite_element* create_sub_element(std::size_t i) const = 0;

    /// Create a new class instance
    virtual finite_element* create() const = 0;

  };

  /// This class defines the interface for a local-to-global mapping of
  /// degrees of freedom (dofs).

  class dofmap
  {
  public:

    /// Destructor
    virtual ~dofmap() {}

    /// Return a string identifying the dofmap
    virtual const char* signature() const = 0;

    /// Return true iff mesh entities of topological dimension d are
    /// needed
    virtual bool needs_mesh_entities(std::size_t d) const = 0;

    /// Return the topological dimension of the associated cell shape
    virtual std::size_t topological_dimension() const = 0;

    /// Return the geometric dimension of the associated cell shape
    virtual std::size_t geometric_dimension() const = 0;

    /// Return the dimension of the global finite element function space
    virtual std::size_t global_dimension(const std::vector<std::size_t>&
                                         num_global_mesh_entities) const = 0;

    /// Return the dimension of the local finite element function space
    /// for a cell
    virtual std::size_t local_dimension() const = 0;

    /// Return the number of dofs on each cell facet
    virtual std::size_t num_facet_dofs() const = 0;

    /// Return the number of dofs associated with each cell entity of
    /// dimension d
    virtual std::size_t num_entity_dofs(std::size_t d) const = 0;

    /// Tabulate the local-to-global mapping of dofs on a cell
    virtual void tabulate_dofs(std::size_t* dofs,
                               const std::vector<std::size_t>& num_global_entities,
                               const cell& c) const = 0;

    /// Tabulate the local-to-local mapping from facet dofs to cell dofs
    virtual void tabulate_facet_dofs(std::size_t* dofs,
                                     std::size_t facet) const = 0;

    /// Tabulate the local-to-local mapping of dofs on entity (d, i)
    virtual void tabulate_entity_dofs(std::size_t* dofs,
                                      std::size_t d, std::size_t i) const = 0;

    /// Tabulate the coordinates of all dofs on a cell
    virtual void tabulate_coordinates(double** dof_coordinates,
                                      const double* vertex_coordinates) const = 0;

    /// Return the number of sub dofmaps (for a mixed element)
    virtual std::size_t num_sub_dofmaps() const = 0;

    /// Create a new dofmap for sub dofmap i (for a mixed element)
    virtual dofmap* create_sub_dofmap(std::size_t i) const = 0;

    /// Create a new class instance
    virtual dofmap* create() const = 0;

  };

  /// This class defines the interface for the tabulation of the cell
  /// tensor corresponding to the local contribution to a form from
  /// the integral over a cell.

  class cell_integral
  {
  public:

    /// Destructor
    virtual ~cell_integral() {}

    /// Tabulate the tensor for the contribution from a local cell
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates,
                                 int cell_orientation) const = 0;

    // FIXME: New experimental version
    /// Tabulate the tensor for the contribution from a local cell
    virtual void tabulate_tensor_new(double* A,
                                     const double * const * w,
                                     const cell_geometry& c) const {}

  };

  /// This class defines the interface for the tabulation of the
  /// exterior facet tensor corresponding to the local contribution to
  /// a form from the integral over an exterior facet.

  class exterior_facet_integral
  {
  public:

    /// Destructor
    virtual ~exterior_facet_integral() {}

    /// Tabulate the tensor for the contribution from a local exterior facet
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates,
                                 std::size_t facet) const = 0;

  };

  /// This class defines the interface for the tabulation of the
  /// interior facet tensor corresponding to the local contribution to
  /// a form from the integral over an interior facet.

  class interior_facet_integral
  {
  public:

    /// Destructor
    virtual ~interior_facet_integral() {}

    /// Tabulate the tensor for the contribution from a local interior facet
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates_0,
                                 const double* vertex_coordinates_1,
                                 std::size_t facet_0,
                                 std::size_t facet_1) const = 0;

  };

  /// This class defines the interface for the tabulation of
  /// an expression evaluated at exactly one point.

  class point_integral
  {
  public:

    /// Constructor
    point_integral() {}

    /// Destructor
    virtual ~point_integral() {}

    /// Tabulate the tensor for the contribution from the local vertex
    virtual void tabulate_tensor(double* A,
                                 const double * const * w,
                                 const double* vertex_coordinates,
                                 std::size_t vertex) const = 0;

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

    /// Destructor
    virtual ~form() {}

    /// Return a string identifying the form
    virtual const char* signature() const = 0;

    /// Return the rank of the global tensor (r)
    virtual std::size_t rank() const = 0;

    /// Return the number of coefficients (n)
    virtual std::size_t num_coefficients() const = 0;

    /// Return the number of cell domains
    virtual std::size_t num_cell_domains() const = 0;

    /// Return the number of exterior facet domains
    virtual std::size_t num_exterior_facet_domains() const = 0;

    /// Return the number of interior facet domains
    virtual std::size_t num_interior_facet_domains() const = 0;

    /// Return the number of point domains
    virtual std::size_t num_point_domains() const = 0;

    /// Return whether form has any cell integrals
    virtual bool has_cell_integrals() const = 0;

    /// Return whether form has any exterior facet integrals
    virtual bool has_exterior_facet_integrals() const = 0;

    /// Return whether form has any interior facet integrals
    virtual bool has_interior_facet_integrals() const = 0;

    /// Return whether form has any point integrals
    virtual bool has_point_integrals() const = 0;

    /// Create a new finite element for argument function i
    virtual finite_element* create_finite_element(std::size_t i) const = 0;

    /// Create a new dofmap for argument function i
    virtual dofmap* create_dofmap(std::size_t i) const = 0;

    /// Create a new cell integral on sub domain i
    virtual cell_integral* create_cell_integral(std::size_t i) const = 0;

    /// Create a new exterior facet integral on sub domain i
    virtual exterior_facet_integral*
    create_exterior_facet_integral(std::size_t i) const = 0;

    /// Create a new interior facet integral on sub domain i
    virtual interior_facet_integral*
    create_interior_facet_integral(std::size_t i) const = 0;

    /// Create a new point integral on sub domain i
    virtual point_integral* create_point_integral(std::size_t i) const = 0;

    /// Create a new cell integral on everywhere else
    virtual cell_integral* create_default_cell_integral() const = 0;

    /// Create a new exterior facet integral on everywhere else
    virtual exterior_facet_integral*
    create_default_exterior_facet_integral() const = 0;

    /// Create a new interior facet integral on everywhere else
    virtual interior_facet_integral*
    create_default_interior_facet_integral() const = 0;

    /// Create a new point integral on everywhere else
    virtual point_integral* create_default_point_integral() const = 0;

  };

}

#endif
