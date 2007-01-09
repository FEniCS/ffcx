"Headers for the UFC 1.0 format."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-08 -- 2007-01-08"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

dof_map = """\
class %s : public ufc::dof_map()
{
public:

  /// Constructor
  %s() : dof_map()
  {
%(hej)s
  }

  /// Destructor
  ~%s()
  {
%s
  }

  /// Return a string identifying the dof map
  const char* signature() const
  {
%s
  }

  /// Return true iff mesh entities of topological dimension d are needed
  bool needs_mesh_entities(unsigned int d) const
  {
%s
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  bool init_mesh(const mesh& mesh)
  {
%s
  }

  /// Initialize dof map for given cell
  void init_cell(const mesh& mesh,
                 const cell& cell)
  {
%s
  }

  /// Return the dimension of the global finite element function space
  unsigned int global_dimension() const
  {
%s
  }

  /// Return the dimension of the local finite element function space
  unsigned int local_dimension() const
  {
%s
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(unsigned int* dofs,
                     const mesh& m,
                     const cell& c) const
  {
%s
  }

  /// Tabulate the local-to-global mapping of dofs on a facet of a cell
  void tabulate_facet_dofs(unsigned int* dofs,
                           const mesh& m,
                           const cell& c,
                           unsigned int facet) const
  {
%s
  }

};
"""

finite_element = """\
class %s : public ufc::finite_element
{
public:

  /// Constructor
  %s : finite_element()
  {
%s
  }

  /// Destructor
  ~%s()
  {
%s
  }

  /// Return a string identifying the finite element
  const char* signature() const
  {
%s
  }

  /// Return the cell shape
  ufc::shape cell_shape() const
  {
%s
  }

  /// Return the dimension of the finite element function space
  unsigned int space_dimension() const
  {
%s
  }

  /// Return the rank of the value space
  unsigned int value_rank() const
  {
%s
  }

  /// Return the dimension of the value space for axis i
  unsigned int value_dimension(unsigned int i) const
  {
%s
  }

  /// Evaluate basis function i at the point x = (x[0], x[1], ...) in cell
  void evaluate_basis(double* values,
                      const double* x,
                      unsigned int i,
                      const ufc::cell& c) const
  {
%s
  }

  /// Evaluate linear functional for dof i on the function f
  double evaluate_dof(unsigned int i,
                      const ufc::function& f,
                      const ufc::cell& c) const
  {
%s
  }

  /// Interpolate vertex values from dof values
  void interpolate_vertex_values(double* vertex_values,
                                 const double* dof_values) const
  {
%s
  }

  /// Return the number of sub elements (for a mixed finite element)
  unsigned int num_sub_elements(unsigned int i) const
  {
%s
  }

  /// Return sub element i (for a mixed finite element)
  const ufc::finite_element& sub_element(unsigned int i) const
  {
%s
  }

};
"""

cell_integral = """\
class %s : public ufc::cell_integral
{
public:

  /// Constructor
  %s : cell_integral()
  {
%s
  }

  /// Destructor
  ~%s()
  {
%s
  }

  /// Tabulate the tensor for the contribution from a local cell
  void tabulate_tensor(double* A,
                       const double * const * w,
                       const cell& c) const
  {
%s
  }

};
"""

exterior_facet_integral = """\
class %s : public ufc::exterior_facet_integral
{
public:

  /// Constructor
  %s : exterior_facet_integral()
  {
%s
  }

  /// Destructor
  ~%s()
  {
%s
  }

  /// Tabulate the tensor for the contribution from a local exterior facet
  void tabulate_tensor(double* A,
                       const double * const * w,
                       const cell& c,
                       unsigned int facet) const
  {
%s
  }

};
"""

interior_facet_integral = """\
class %s : public ufc::interior_facet_integral
{
public:

  /// Constructor
  %s : interior_facet_integral()
  {
%s
  }

  /// Destructor
  ~%s()
  {
%s
  }

  /// Tabulate the tensor for the contribution from a local interior facet
  void tabulate_tensor(double* A,
                       const double * const * w,
                       const cell& c0,
                       const cell& c1,
                       unsigned int facet0,
                       unsigned int facet1) const
  {
%s
  }

};
"""

form = """\
class %s : public ufc::form
{
public:

  /// Constructor
  %s : form()
  {
%s
  }

  /// Destructor
  ~%s()
  {
%s
  }

  /// Return a string identifying the form
  const char* signature() const
  {
%s
  }

  /// Return the rank of the element tensor (r)
  unsigned int rank() const
  {
%s
  }

  /// Return the number of coefficients (n)
  unsigned int num_coefficients() const
  {
%s
  }

  /// Create a new dof map for argument function i
  dof_map* create_dof_map(unsigned int i) const
  {
%s
  }

  /// Create a new finite element for argument function i
  finite_element* create_finite_element(unsigned int i) const
  {
%s
  }

  /// Create a new cell integral (return 0 if contribution is zero)
  cell_integral* create_cell_integral() const
  {
%s
  }

  /// Create a new interior facet integral (return 0 if contribution is zero)
  interior_facet_integral* create_interior_facet_integral() const
  {
%s
  }

  /// Create a new exterior facet integral (return 0 if contribution is zero)
  exterior_facet_integral* create_exterior_facet_integral() const
  {
%s
  }
  
};
"""
