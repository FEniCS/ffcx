# Code generation format strings for UFC (Unified Form-assembly Code) v. 2.3.0+.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2014

cell_integral_combined = """\
/// This class defines the interface for the tabulation of the cell
/// tensor corresponding to the local contribution to a form from
/// the integral over a cell.

class %(classname)s: public ufc::cell_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::cell_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from a local cell
  virtual void tabulate_tensor(double* %(restrict)s A,
                               const double * const * %(restrict)s w,
                               const double* %(restrict)s vertex_coordinates,
                               int cell_orientation) const
  {
%(tabulate_tensor)s
  }

};
"""

cell_integral_header = """\
/// This class defines the interface for the tabulation of the cell
/// tensor corresponding to the local contribution to a form from
/// the integral over a cell.

class %(classname)s: public ufc::cell_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  virtual ~%(classname)s();

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const;

  /// Tabulate the tensor for the contribution from a local cell
  virtual void tabulate_tensor(double* %(restrict)s A,
                               const double * const * %(restrict)s w,
                               const double* %(restrict)s vertex_coordinates,
                               int cell_orientation) const;

};
"""

cell_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::cell_integral()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Tabulate which form coefficients are used by this integral
const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

/// Tabulate the tensor for the contribution from a local cell
void %(classname)s::tabulate_tensor(double* %(restrict)s A,
                                    const double * const * %(restrict)s w,
                                    const double* %(restrict)s vertex_coordinates,
                                    int cell_orientation) const
{
%(tabulate_tensor)s
}
"""

exterior_facet_integral_combined = """\
/// This class defines the interface for the tabulation of the
/// exterior facet tensor corresponding to the local contribution to
/// a form from the integral over an exterior facet.

class %(classname)s: public ufc::exterior_facet_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::exterior_facet_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from a local exterior facet
  virtual void tabulate_tensor(double* %(restrict)s A,
                               const double * const * %(restrict)s w,
                               const double* %(restrict)s vertex_coordinates,
                               std::size_t facet,
                               int cell_orientation) const
  {
%(tabulate_tensor)s
  }

};
"""

exterior_facet_integral_header = """\
/// This class defines the interface for the tabulation of the
/// exterior facet tensor corresponding to the local contribution to
/// a form from the integral over an exterior facet.

class %(classname)s: public ufc::exterior_facet_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  virtual ~%(classname)s();

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const;

  /// Tabulate the tensor for the contribution from a local exterior facet
  virtual void tabulate_tensor(double* %(restrict)s A,
                               const double * const * %(restrict)s w,
                               const double* %(restrict)s vertex_coordinates,
                               std::size_t facet,
                               int cell_orientation) const;

};
"""

exterior_facet_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::exterior_facet_integral()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Tabulate which form coefficients are used by this integral
const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

/// Tabulate the tensor for the contribution from a local exterior facet
void %(classname)s::tabulate_tensor(double* %(restrict)s A,
                                    const double * const * %(restrict)s w,
                                    const double* %(restrict)s vertex_coordinates,
                                    std::size_t facet,
                                    int cell_orientation) const
{
%(tabulate_tensor)s
}
"""

interior_facet_integral_combined = """\
/// This class defines the interface for the tabulation of the
/// interior facet tensor corresponding to the local contribution to
/// a form from the integral over an interior facet.

class %(classname)s: public ufc::interior_facet_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::interior_facet_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from a local interior facet
  virtual void tabulate_tensor(double* %(restrict)s A,
                               const double * const * %(restrict)s w,
                               const double* %(restrict)s vertex_coordinates_0,
                               const double* %(restrict)s vertex_coordinates_1,
                               std::size_t facet_0,
                               std::size_t facet_1,
                               int cell_orientation_0,
                               int cell_orientation_1) const
  {
%(tabulate_tensor)s
  }

};
"""

interior_facet_integral_header = """\
/// This class defines the interface for the tabulation of the
/// interior facet tensor corresponding to the local contribution to
/// a form from the integral over an interior facet.

class %(classname)s: public ufc::interior_facet_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  virtual ~%(classname)s();

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const;

  /// Tabulate the tensor for the contribution from a local interior facet
  virtual void tabulate_tensor(double* %(restrict)s A,
                               const double * const * %(restrict)s w,
                               const double* %(restrict)s vertex_coordinates_0,
                               const double* %(restrict)s vertex_coordinates_1,
                               std::size_t facet_0,
                               std::size_t facet_1,
                               int cell_orientation_0,
                               int cell_orientation_1) const;

};
"""

interior_facet_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::interior_facet_integral()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Tabulate which form coefficients are used by this integral
const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

/// Tabulate the tensor for the contribution from a local interior facet
void %(classname)s::tabulate_tensor(double* %(restrict)s A,
                                    const double * const * %(restrict)s w,
                                    const double* %(restrict)s vertex_coordinates_0,
                                    const double* %(restrict)s vertex_coordinates_1,
                                    std::size_t facet_0,
                                    std::size_t facet_1,
                                    int cell_orientation_0,
                                    int cell_orientation_1) const
{
%(tabulate_tensor)s
}
"""

point_integral_combined = """\
/// This class defines the interface for the tabulation of
/// an expression evaluated at exactly one point.

class %(classname)s: public ufc::point_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::point_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from the local vertex
  virtual void tabulate_tensor(double* %(restrict)s A,
                               const double * const * %(restrict)s w,
                               const double* %(restrict)s vertex_coordinates,
                               std::size_t vertex,
                               int cell_orientation) const
  {
%(tabulate_tensor)s
  }

};
"""

point_integral_header = """\
/// This class defines the interface for the tabulation of
/// an expression evaluated at exactly one point.

class %(classname)s: public ufc::point_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  virtual ~%(classname)s();

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const;

  /// Tabulate the tensor for the contribution from the local vertex
  virtual void tabulate_tensor(double* %(restrict)s A,
                               const double * const * %(restrict)s w,
                               const double* %(restrict)s vertex_coordinates,
                               std::size_t vertex,
                               int cell_orientation) const;

};
"""

point_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::point_integral()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Tabulate which form coefficients are used by this integral
const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

/// Tabulate the tensor for the contribution from the local vertex
void %(classname)s::tabulate_tensor(double* %(restrict)s A,
                                    const double * const * %(restrict)s w,
                                    const double* %(restrict)s vertex_coordinates,
                                    std::size_t vertex,
                                    int cell_orientation) const
{
%(tabulate_tensor)s
}
"""

custom_integral_combined = """\
/// This class defines the interface for the tabulation of the
/// tensor corresponding to the local contribution to a form from
/// the integral over a custom domain defined in terms of a set of
/// quadrature points and weights.

class %(classname)s: public ufc::custom_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::custom_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const
  {
%(enabled_coefficients)s
  }

  /// Return the number of cells involved in evaluation of the integral
  virtual std::size_t num_cells() const
  {
%(num_cells)s
  }

  /// Tabulate the tensor for the contribution from custom domain
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const double* vertex_coordinates,
                               std::size_t num_quadrature_points,
                               const double* quadrature_points,
                               const double* quadrature_weights,
                               int cell_orientation) const
  {
%(tabulate_tensor)s
  }

};
"""

custom_integral_header = """\
/// This class defines the interface for the tabulation of the
/// tensor corresponding to the local contribution to a form from
/// the integral over a custom domain defined in terms of a set of
/// quadrature points and weights.

class %(classname)s: public ufc::custom_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  virtual ~%(classname)s();

  /// Tabulate which form coefficients are used by this integral
  virtual const std::vector<bool> & enabled_coefficients() const;

  /// Return the number of cells involved in evaluation of the integral
  virtual std::size_t num_cells() const;

  /// Tabulate the tensor for the contribution from custom domain
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const double* vertex_coordinates,
                               std::size_t num_quadrature_points,
                               const double* quadrature_points,
                               const double* quadrature_weights,
                               int cell_orientation) const;

};
"""

custom_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::custom_integral()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Tabulate which form coefficients are used by this integral
const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

/// Return the number of cells involved in evaluation of the integral
std::size_t %(classname)s::num_cells() const
{
%(num_cells)s
}

/// Tabulate the tensor for the contribution from custom domain
void %(classname)s::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const double* vertex_coordinates,
                                    std::size_t num_quadrature_points,
                                    const double* quadrature_points,
                                    const double* quadrature_weights,
                                    int cell_orientation) const
{
%(tabulate_tensor)s
}
"""
