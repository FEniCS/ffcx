# Code generation format strings for UFC (Unified Form-assembly Code) v. 2.0.5.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2011.

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
  virtual unsigned int topological_dimension() const
  {
%(topological_dimension)s
  }

  /// Return the geometric dimension of the cell shape
  virtual unsigned int geometric_dimension() const
  {
%(geometric_dimension)s
  }

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const
  {
%(space_dimension)s
  }

  /// Return the rank of the value space
  virtual unsigned int value_rank() const
  {
%(value_rank)s
  }

  /// Return the dimension of the value space for axis i
  virtual unsigned int value_dimension(unsigned int i) const
  {
%(value_dimension)s
  }

  /// Evaluate basis function i at given point in cell
  virtual void evaluate_basis(unsigned int i,
                              double* values,
                              const double* coordinates,
                              const ufc::cell& c) const
  {
%(evaluate_basis)s
  }

  /// Evaluate all basis functions at given point in cell
  virtual void evaluate_basis_all(double* values,
                                  const double* coordinates,
                                  const ufc::cell& c) const
  {
%(evaluate_basis_all)s
  }

  /// Evaluate order n derivatives of basis function i at given point in cell
  virtual void evaluate_basis_derivatives(unsigned int i,
                                          unsigned int n,
                                          double* values,
                                          const double* coordinates,
                                          const ufc::cell& c) const
  {
%(evaluate_basis_derivatives)s
  }

  /// Evaluate order n derivatives of all basis functions at given point in cell
  virtual void evaluate_basis_derivatives_all(unsigned int n,
                                              double* values,
                                              const double* coordinates,
                                              const ufc::cell& c) const
  {
%(evaluate_basis_derivatives_all)s
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
%(evaluate_dof)s
  }

  /// Evaluate linear functionals for all dofs on the function f
  virtual void evaluate_dofs(double* values,
                             const ufc::function& f,
                             const ufc::cell& c) const
  {
%(evaluate_dofs)s
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const
  {
%(interpolate_vertex_values)s
  }

  /// Map coordinate xhat from reference cell to coordinate x in cell
  virtual void map_from_reference_cell(double* x,
                                       const double* xhat,
                                       const ufc::cell& c) const
  {
%(map_from_reference_cell)s
  }

  /// Map from coordinate x in cell to coordinate xhat in reference cell
  virtual void map_to_reference_cell(double* xhat,
                                     const double* x,
                                     const ufc::cell& c) const
  {
%(map_to_reference_cell)s
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
%(num_sub_elements)s
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
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
  virtual unsigned int topological_dimension() const;

  /// Return the geometric dimension of the cell shape
  virtual unsigned int geometric_dimension() const;

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const;

  /// Return the rank of the value space
  virtual unsigned int value_rank() const;

  /// Return the dimension of the value space for axis i
  virtual unsigned int value_dimension(unsigned int i) const;

  /// Evaluate basis function i at given point in cell
  virtual void evaluate_basis(unsigned int i,
                              double* values,
                              const double* coordinates,
                              const ufc::cell& c) const;

  /// Evaluate all basis functions at given point in cell
  virtual void evaluate_basis_all(double* values,
                                  const double* coordinates,
                                  const ufc::cell& c) const;

  /// Evaluate order n derivatives of basis function i at given point in cell
  virtual void evaluate_basis_derivatives(unsigned int i,
                                          unsigned int n,
                                          double* values,
                                          const double* coordinates,
                                          const ufc::cell& c) const;

  /// Evaluate order n derivatives of all basis functions at given point in cell
  virtual void evaluate_basis_derivatives_all(unsigned int n,
                                              double* values,
                                              const double* coordinates,
                                              const ufc::cell& c) const;

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const;

  /// Evaluate linear functionals for all dofs on the function f
  virtual void evaluate_dofs(double* values,
                             const ufc::function& f,
                             const ufc::cell& c) const;

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const;

  /// Map coordinate xhat from reference cell to coordinate x in cell
  virtual void map_from_reference_cell(double* x,
                                       const double* xhat,
                                       const ufc::cell& c) const;

  /// Map from coordinate x in cell to coordinate xhat in reference cell
  virtual void map_to_reference_cell(double* xhat,
                                     const double* x,
                                     const ufc::cell& c) const;

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const;

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const;

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
unsigned int %(classname)s::topological_dimension() const
{
%(topological_dimension)s
}

/// Return the geometric dimension of the cell shape
unsigned int %(classname)s::geometric_dimension() const
{
%(geometric_dimension)s
}

/// Return the dimension of the finite element function space
unsigned int %(classname)s::space_dimension() const
{
%(space_dimension)s
}

/// Return the rank of the value space
unsigned int %(classname)s::value_rank() const
{
%(value_rank)s
}

/// Return the dimension of the value space for axis i
unsigned int %(classname)s::value_dimension(unsigned int i) const
{
%(value_dimension)s
}

/// Evaluate basis function i at given point in cell
void %(classname)s::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
%(evaluate_basis)s
}

/// Evaluate all basis functions at given point in cell
void %(classname)s::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
%(evaluate_basis_all)s
}

/// Evaluate order n derivatives of basis function i at given point in cell
void %(classname)s::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
%(evaluate_basis_derivatives)s
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void %(classname)s::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
%(evaluate_basis_derivatives_all)s
}

/// Evaluate linear functional for dof i on the function f
double %(classname)s::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
%(evaluate_dof)s
}

/// Evaluate linear functionals for all dofs on the function f
void %(classname)s::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
%(evaluate_dofs)s
}

/// Interpolate vertex values from dof values
void %(classname)s::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
%(interpolate_vertex_values)s
}

/// Map coordinate xhat from reference cell to coordinate x in cell
void %(classname)s::map_from_reference_cell(double* x,
                                            const double* xhat,
                                            const ufc::cell& c) const
{
%(map_from_reference_cell)s
}

/// Map from coordinate x in cell to coordinate xhat in reference cell
void %(classname)s::map_to_reference_cell(double* xhat,
                                          const double* x,
                                          const ufc::cell& c) const
{
%(map_to_reference_cell)s
}

/// Return the number of sub elements (for a mixed element)
unsigned int %(classname)s::num_sub_elements() const
{
%(num_sub_elements)s
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* %(classname)s::create_sub_element(unsigned int i) const
{
%(create_sub_element)s
}

/// Create a new class instance
ufc::finite_element* %(classname)s::create() const
{
%(create)s
}
"""
