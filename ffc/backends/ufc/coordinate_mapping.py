# Code generation format strings for UFC (Unified Form-assembly Code) v. 1.7.0dev.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2015.

coordinate_mapping_header = """\
/// A representation of a coordinate mapping parameterized by a local finite element basis on each cell
class %(classname)s: public ufc::coordinate_mapping
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  ~%(classname)s() override;

  /// Return coordinate_mapping signature string
  const char * signature() const final override;

  /// Create object of the same type
  coordinate_mapping * create() const final override;

  /// Return geometric dimension of the coordinate_mapping
  std::size_t geometric_dimension() const final override;

  /// Return topological dimension of the coordinate_mapping
  std::size_t topological_dimension() const final override;

  /// Return cell shape of the coordinate_mapping
  ufc::shape cell_shape() const final override;

  /// Create finite_element object representing the coordinate parameterization
  ufc::finite_element * create_coordinate_finite_element() const final override;

  /// Create dofmap object representing the coordinate parameterization
  ufc::dofmap * create_coordinate_dofmap() const final override;

  /// Compute physical coordinates x from reference coordinates X, the inverse of compute_reference_coordinates
  void compute_physical_coordinates(
      double * x, std::size_t num_points,
      const double * X,
      const double * coordinate_dofs, int cell_orientation) const final override;

  /// Compute reference coordinates X from physical coordinates x, the inverse of compute_physical_coordinates
  void compute_reference_coordinates(
      double * X, std::size_t num_points,
      const double * x,
      const double * coordinate_dofs, int cell_orientation) const final override;

  /// Compute Jacobian of coordinate mapping J = dx/dX at reference coordinates X
  void compute_jacobians(
      double * J, std::size_t num_points,
      const double * X,
      const double * coordinate_dofs, int cell_orientation) const final override;

  /// Compute determinants of (pseudo-)Jacobians J
  void compute_jacobian_determinants(
      double * detJ, std::size_t num_points,
      const double * J) const final override;

  /// Compute (pseudo-)inverses K of (pseudo-)Jacobians J
  void compute_jacobian_inverses(
      double * K, std::size_t num_points,
      const double * J, const double * detJ) const final override;

  /// Combined (for convenience) computation of x, J, detJ, K from X and coordinate_dofs on a cell
  void compute_geometry(
      double * x, double * J, double * detJ, double * K, std::size_t num_points,
      const double * X,
      const double * coordinate_dofs, int cell_orientation) const final override;
};
"""

coordinate_mapping_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::coordinate_mapping()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s() override
{
%(destructor)s
}

/// Return coordinate_mapping signature string
const char * %(classname)s::signature() const final override
{
%(signature)s
}

/// Create object of the same type
coordinate_mapping * %(classname)s::create() const final override
{
%(create)s
}

/// Return geometric dimension of the coordinate_mapping
std::size_t %(classname)s::geometric_dimension() const final override
{
%(geometric_dimension)s
}

/// Return topological dimension of the coordinate_mapping
std::size_t %(classname)s::topological_dimension() const final override
{
%(topological_dimension)s
}

/// Return cell shape of the coordinate_mapping
ufc::shape %(classname)s::cell_shape() const final override
{
%(cell_shape)s
}

/// Create finite_element object representing the coordinate parameterization
ufc::finite_element * %(classname)s::create_coordinate_finite_element() const final override
{
%(create_coordinate_finite_element)s
}

/// Create dofmap object representing the coordinate parameterization
ufc::dofmap * %(classname)s::create_coordinate_dofmap() const final override
{
%(create_coordinate_dofmap)s
}

/// Compute physical coordinates x from reference coordinates X, the inverse of compute_reference_coordinates
void %(classname)s::compute_physical_coordinates(
    double * x, std::size_t num_points,
    const double * X,
    const double * coordinate_dofs, int cell_orientation) const final override
{
%(compute_physical_coordinates)s
}

/// Compute reference coordinates X from physical coordinates x, the inverse of compute_physical_coordinates
void %(classname)s::compute_reference_coordinates(
    double * X, std::size_t num_points,
    const double * x,
    const double * coordinate_dofs, int cell_orientation) const final override
{
%(compute_reference_coordinates)s
}

/// Compute Jacobian of coordinate mapping J = dx/dX at reference coordinates X
void %(classname)s::compute_jacobians(
    double * J, std::size_t num_points,
    const double * X,
    const double * coordinate_dofs, int cell_orientation) const final override
{
%(compute_jacobians)s
}

/// Compute determinants of (pseudo-)Jacobians J
void %(classname)s::compute_jacobian_determinants(
    double * detJ, std::size_t num_points,
    const double * J) const final override
{
%(compute_jacobian_determinants)s
}

/// Compute (pseudo-)inverses K of (pseudo-)Jacobians J
void %(classname)s::compute_jacobian_inverses(
    double * K, std::size_t num_points,
    const double * J, const double * detJ) const final override
{
%(compute_jacobian_inverses)s
}

/// Combined (for convenience) computation of x, J, detJ, K from X and coordinate_dofs on a cell
void %(classname)s::compute_geometry(
    double * x, double * J, double * detJ, double * K, std::size_t num_points,
    const double * X,
    const double * coordinate_dofs, int cell_orientation) const final override
{
%(compute_geometry)s
}
"""

coordinate_mapping_combined = """\
/// A representation of a geometric domain parameterized by a local basis on each cell
class %(classname)s: public ufc::coordinate_mapping
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::coordinate_mapping()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Return coordinate_mapping signature string
  const char * signature() const final override
  {
%(signature)s
  }

  /// Create object of the same type
  coordinate_mapping * create() const final override
  {
%(create)s
  }

  /// Return geometric dimension of the coordinate_mapping
  std::size_t geometric_dimension() const final override
  {
%(geometric_dimension)s
  }

  /// Return topological dimension of the coordinate_mapping
  std::size_t topological_dimension() const final override
  {
%(topological_dimension)s
  }

  /// Return cell shape of the coordinate_mapping
  ufc::shape cell_shape() const final override
  {
%(cell_shape)s
  }

  /// Create finite_element object representing the coordinate parameterization
  ufc::finite_element * create_coordinate_finite_element() const final override
  {
%(create_coordinate_finite_element)s
  }

  /// Create dofmap object representing the coordinate parameterization
  ufc::dofmap * create_coordinate_dofmap() const final override
  {
%(create_coordinate_dofmap)s
  }

  /// Compute physical coordinates x from reference coordinates X, the inverse of compute_reference_coordinates
  void compute_physical_coordinates(
      double * x, std::size_t num_points,
      const double * X,
      const double * coordinate_dofs, int cell_orientation) const final override
  {
%(compute_physical_coordinates)s
  }

  /// Compute reference coordinates X from physical coordinates x, the inverse of compute_physical_coordinates
  void compute_reference_coordinates(
      double * X, std::size_t num_points,
      const double * x,
      const double * coordinate_dofs, int cell_orientation) const final override
  {
%(compute_reference_coordinates)s
  }

  /// Compute Jacobian of coordinate mapping J = dx/dX at reference coordinates X
  void compute_jacobians(
      double * J, std::size_t num_points,
      const double * X,
      const double * coordinate_dofs, int cell_orientation) const final override
  {
%(compute_jacobians)s
  }

  /// Compute determinants of (pseudo-)Jacobians J
  void compute_jacobian_determinants(
      double * detJ, std::size_t num_points,
      const double * J) const final override
  {
%(compute_jacobian_determinants)s
  }

  /// Compute (pseudo-)inverses K of (pseudo-)Jacobians J
  void compute_jacobian_inverses(
      double * K, std::size_t num_points,
      const double * J, const double * detJ) const final override
  {
%(compute_jacobian_inverses)s
  }

  /// Combined (for convenience) computation of x, J, detJ, K from X and coordinate_dofs on a cell
  void compute_geometry(
      double * x, double * J, double * detJ, double * K, std::size_t num_points,
      const double * X,
      const double * coordinate_dofs, int cell_orientation) const final override
  {
%(compute_geometry)s
  }

};
"""
