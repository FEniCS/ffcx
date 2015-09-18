# Code generation format strings for UFC (Unified Form-assembly Code) v. 1.7.0dev.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2015.

finite_element_combined = """\
/// This class defines the interface for a finite element.

class %(classname)s: public ufc::finite_element
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::finite_element()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
%(signature)s
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
%(cell_shape)s
  }

  /// Return the topological dimension of the cell shape
  virtual std::size_t topological_dimension() const
  {
%(topological_dimension)s
  }

  /// Return the geometric dimension of the cell shape
  virtual std::size_t geometric_dimension() const
  {
%(geometric_dimension)s
  }

  /// Return the dimension of the finite element function space
  virtual std::size_t space_dimension() const
  {
%(space_dimension)s
  }

  /// Return the rank of the value space
  virtual std::size_t value_rank() const
  {
%(value_rank)s
  }

  /// Return the dimension of the value space for axis i
  virtual std::size_t value_dimension(std::size_t i) const
  {
%(value_dimension)s
  }

  /// Return the number of components of the value space
  virtual std::size_t value_size() const
  {
%(value_size)s
  }

  /// Return the rank of the reference value space
  virtual std::size_t reference_value_rank() const
  {
%(reference_value_rank)s
  }

  /// Return the dimension of the reference value space for axis i
  virtual std::size_t reference_value_dimension(std::size_t i) const
  {
%(reference_value_dimension)s
  }

  /// Return the number of components of the reference value space
  virtual std::size_t reference_value_size() const
  {
%(reference_value_size)s
  }

  /// Evaluate basis function i at given point x in cell (actual implementation)
  static void _evaluate_basis(std::size_t i,
                              double* values,
                              const double* x,
                              const double* coordinate_dofs,
                              int cell_orientation)
  {
%(evaluate_basis)s
  }

  /// Evaluate basis function i at given point x in cell (non-static member function)
  virtual void evaluate_basis(std::size_t i,
                              double* values,
                              const double* x,
                              const double* coordinate_dofs,
                              int cell_orientation) const
  {
    _evaluate_basis(i, values, x, coordinate_dofs, cell_orientation);
  }

  /// Evaluate all basis functions at given point x in cell (actual implementation)
  static void _evaluate_basis_all(double* values,
                                  const double* x,
                                  const double* coordinate_dofs,
                                  int cell_orientation)
  {
%(evaluate_basis_all)s
  }

  /// Evaluate all basis functions at given point x in cell (non-static member function)
  virtual void evaluate_basis_all(double* values,
                                  const double* x,
                                  const double* coordinate_dofs,
                                  int cell_orientation) const
  {
    _evaluate_basis_all(values, x, coordinate_dofs, cell_orientation);
  }

  /// Evaluate order n derivatives of basis function i at given point x in cell (actual implementation)
  static void _evaluate_basis_derivatives(std::size_t i,
                                          std::size_t n,
                                          double* values,
                                          const double* x,
                                          const double* coordinate_dofs,
                                          int cell_orientation)
  {
%(evaluate_basis_derivatives)s
  }

  /// Evaluate order n derivatives of basis function i at given point x in cell (non-static member function)
  virtual void evaluate_basis_derivatives(std::size_t i,
                                          std::size_t n,
                                          double* values,
                                          const double* x,
                                          const double* coordinate_dofs,
                                          int cell_orientation) const
  {
    _evaluate_basis_derivatives(i, n, values, x, coordinate_dofs, cell_orientation);
  }

  /// Evaluate order n derivatives of all basis functions at given point x in cell (actual implementation)
  static void _evaluate_basis_derivatives_all(std::size_t n,
                                              double* values,
                                              const double* x,
                                              const double* coordinate_dofs,
                                              int cell_orientation)
  {
%(evaluate_basis_derivatives_all)s
  }

  /// Evaluate order n derivatives of all basis functions at given point x in cell (non-static member function)
  virtual void evaluate_basis_derivatives_all(std::size_t n,
                                              double* values,
                                              const double* x,
                                              const double* coordinate_dofs,
                                              int cell_orientation) const
  {
    _evaluate_basis_derivatives_all(n, values, x, coordinate_dofs, cell_orientation);
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(std::size_t i,
                              const ufc::function& f,
                              const double* coordinate_dofs,
                              int cell_orientation,
                              const ufc::cell& c) const
  {
%(evaluate_dof)s
  }

  /// Evaluate linear functionals for all dofs on the function f
  virtual void evaluate_dofs(double* values,
                             const ufc::function& f,
                             const double* coordinate_dofs,
                             int cell_orientation,
                             const ufc::cell& c) const
  {
%(evaluate_dofs)s
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const double* coordinate_dofs,
                                         int cell_orientation,
                                         const ufc::cell& c) const
  {
%(interpolate_vertex_values)s
  }

  /// Return the number of sub elements (for a mixed element)
  virtual std::size_t num_sub_elements() const
  {
%(num_sub_elements)s
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(std::size_t i) const
  {
%(create_sub_element)s
  }

  /// Create a new class instance
  virtual ufc::finite_element* create() const
  {
%(create)s
  }

};
"""

finite_element_header = """\
/// This class defines the interface for a finite element.

class %(classname)s: public ufc::finite_element
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  virtual ~%(classname)s();

  /// Return a string identifying the finite element
  virtual const char* signature() const;

  /// Return the cell shape
  virtual ufc::shape cell_shape() const;

  /// Return the topological dimension of the cell shape
  virtual std::size_t topological_dimension() const;

  /// Return the geometric dimension of the cell shape
  virtual std::size_t geometric_dimension() const;

  /// Return the dimension of the finite element function space
  virtual std::size_t space_dimension() const;

  /// Return the rank of the value space
  virtual std::size_t value_rank() const;

  /// Return the dimension of the value space for axis i
  virtual std::size_t value_dimension(std::size_t i) const;

  /// Return the number of components of the value space
  virtual std::size_t value_size() const;

  /// Return the rank of the reference value space
  virtual std::size_t reference_value_rank() const;

  /// Return the dimension of the reference value space for axis i
  virtual std::size_t reference_value_dimension(std::size_t i) const;

  /// Return the number of components of the reference value space
  virtual std::size_t reference_value_size() const;

  /// Evaluate basis function i at given point x in cell (actual implementation)
  static void _evaluate_basis(std::size_t i,
                              double* values,
                              const double* x,
                              const double* coordinate_dofs,
                              int cell_orientation);

  /// Evaluate basis function i at given point x in cell (non-static member function)
  virtual void evaluate_basis(std::size_t i,
                              double* values,
                              const double* x,
                              const double* coordinate_dofs,
                              int cell_orientation) const
  {
    _evaluate_basis(i, values, x, coordinate_dofs, cell_orientation);
  }

  /// Evaluate all basis functions at given point x in cell (actual implementation)
  static void _evaluate_basis_all(double* values,
                                  const double* x,
                                  const double* coordinate_dofs,
                                  int cell_orientation);

  /// Evaluate all basis functions at given point x in cell (non-static member function)
  virtual void evaluate_basis_all(double* values,
                                  const double* x,
                                  const double* coordinate_dofs,
                                  int cell_orientation) const
  {
    _evaluate_basis_all(values, x, coordinate_dofs, cell_orientation);
  }

  /// Evaluate order n derivatives of basis function i at given point x in cell (actual implementation)
  static void _evaluate_basis_derivatives(std::size_t i,
                                          std::size_t n,
                                          double* values,
                                          const double* x,
                                          const double* coordinate_dofs,
                                          int cell_orientation);

  /// Evaluate order n derivatives of basis function i at given point x in cell (non-static member function)
  virtual void evaluate_basis_derivatives(std::size_t i,
                                          std::size_t n,
                                          double* values,
                                          const double* x,
                                          const double* coordinate_dofs,
                                          int cell_orientation) const
  {
    _evaluate_basis_derivatives(i, n, values, x, coordinate_dofs, cell_orientation);
  }

  /// Evaluate order n derivatives of all basis functions at given point x in cell (actual implementation)
  static void _evaluate_basis_derivatives_all(std::size_t n,
                                              double* values,
                                              const double* x,
                                              const double* coordinate_dofs,
                                              int cell_orientation);

  /// Evaluate order n derivatives of all basis functions at given point x in cell (non-static member function)
  virtual void evaluate_basis_derivatives_all(std::size_t n,
                                              double* values,
                                              const double* x,
                                              const double* coordinate_dofs,
                                              int cell_orientation) const
  {
    _evaluate_basis_derivatives_all(n, values, x, coordinate_dofs, cell_orientation);
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(std::size_t i,
                              const ufc::function& f,
                              const double* coordinate_dofs,
                              int cell_orientation,
                              const ufc::cell& c) const;

  /// Evaluate linear functionals for all dofs on the function f
  virtual void evaluate_dofs(double* values,
                             const ufc::function& f,
                             const double* coordinate_dofs,
                             int cell_orientation,
                             const ufc::cell& c) const;

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const double* coordinate_dofs,
                                         int cell_orientation,
                                         const ufc::cell& c) const;

  /// Return the number of sub elements (for a mixed element)
  virtual std::size_t num_sub_elements() const;

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(std::size_t i) const;

  /// Create a new class instance
  virtual ufc::finite_element* create() const;

};
"""

finite_element_implementation= """\

/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::finite_element()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Return a string identifying the finite element
const char* %(classname)s::signature() const
{
%(signature)s
}

/// Return the cell shape
ufc::shape %(classname)s::cell_shape() const
{
%(cell_shape)s
}

/// Return the topological dimension of the cell shape
std::size_t %(classname)s::topological_dimension() const
{
%(topological_dimension)s
}

/// Return the geometric dimension of the cell shape
std::size_t %(classname)s::geometric_dimension() const
{
%(geometric_dimension)s
}

/// Return the dimension of the finite element function space
std::size_t %(classname)s::space_dimension() const
{
%(space_dimension)s
}

/// Return the rank of the value space
std::size_t %(classname)s::value_rank() const
{
%(value_rank)s
}

/// Return the dimension of the value space for axis i
std::size_t %(classname)s::value_dimension(std::size_t i) const
{
%(value_dimension)s
}

/// Return the number of components of the value space
virtual std::size_t %(classname)s::value_size() const
{
%(value_size)s
}

/// Return the rank of the reference value space
virtual std::size_t %(classname)s::reference_value_rank() const
{
%(reference_value_rank)s
}

/// Return the dimension of the reference value space for axis i
virtual std::size_t %(classname)s::reference_value_dimension(std::size_t i) const
{
%(reference_value_dimension)s
}

/// Return the number of components of the reference value space
virtual std::size_t %(classname)s::reference_value_size() const
{
%(reference_value_size)s
}

/// Evaluate basis function i at given point x in cell
void %(classname)s::_evaluate_basis(std::size_t i,
                                    double* values,
                                    const double* x,
                                    const double* coordinate_dofs,
                                    int cell_orientation)
{
%(evaluate_basis)s
}

/// Evaluate all basis functions at given point x in cell
void %(classname)s::_evaluate_basis_all(double* values,
                                        const double* x,
                                        const double* coordinate_dofs,
                                        int cell_orientation)
{
%(evaluate_basis_all)s
}

/// Evaluate order n derivatives of basis function i at given point x in cell
void %(classname)s::_evaluate_basis_derivatives(std::size_t i,
                                                std::size_t n,
                                                double* values,
                                                const double* x,
                                                const double* coordinate_dofs,
                                                int cell_orientation)
{
%(evaluate_basis_derivatives)s
}

/// Evaluate order n derivatives of all basis functions at given point x in cell
void %(classname)s::_evaluate_basis_derivatives_all(std::size_t n,
                                                    double* values,
                                                    const double* x,
                                                    const double* coordinate_dofs,
                                                    int cell_orientation)
{
%(evaluate_basis_derivatives_all)s
}

/// Evaluate linear functional for dof i on the function f
double %(classname)s::evaluate_dof(std::size_t i,
                                   const ufc::function& f,
                                   const double* coordinate_dofs,
                                   int cell_orientation,
                                   const ufc::cell& c) const
{
%(evaluate_dof)s
}

/// Evaluate linear functionals for all dofs on the function f
void %(classname)s::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const double* coordinate_dofs,
                                  int cell_orientation,
                                  const ufc::cell& c) const
{
%(evaluate_dofs)s
}

/// Interpolate vertex values from dof values
void %(classname)s::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const double* coordinate_dofs,
                                              int cell_orientation,
                                              const ufc::cell& c) const
{
%(interpolate_vertex_values)s
}

/// Return the number of sub elements (for a mixed element)
std::size_t %(classname)s::num_sub_elements() const
{
%(num_sub_elements)s
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* %(classname)s::create_sub_element(std::size_t i) const
{
%(create_sub_element)s
}

/// Create a new class instance
ufc::finite_element* %(classname)s::create() const
{
%(create)s
}
"""
