# Code generation format strings for UFC (Unified Form-assembly Code) v. 1.7.0dev.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2015

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
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from a local cell
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       int cell_orientation) const final override
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
  ~%(classname)s() override;

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override;

  /// Tabulate the tensor for the contribution from a local cell
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       int cell_orientation) const final override;

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
void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
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
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from a local exterior facet
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t facet,
                       int cell_orientation) const final override
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
  ~%(classname)s() override;

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override;

  /// Tabulate the tensor for the contribution from a local exterior facet
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t facet,
                       int cell_orientation) const final override;

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
void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
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
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from a local interior facet
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs_0,
                       const double * coordinate_dofs_1,
                       std::size_t facet_0,
                       std::size_t facet_1,
                       int cell_orientation_0,
                       int cell_orientation_1) const final override
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
  ~%(classname)s() override;

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override;

  /// Tabulate the tensor for the contribution from a local interior facet
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs_0,
                       const double * coordinate_dofs_1,
                       std::size_t facet_0,
                       std::size_t facet_1,
                       int cell_orientation_0,
                       int cell_orientation_1) const final override;

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
void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs_0,
                                    const double * coordinate_dofs_1,
                                    std::size_t facet_0,
                                    std::size_t facet_1,
                                    int cell_orientation_0,
                                    int cell_orientation_1) const
{
%(tabulate_tensor)s
}
"""

vertex_integral_combined = """\
/// This class defines the interface for the tabulation of
/// an expression evaluated at exactly one vertex.

class %(classname)s: public ufc::vertex_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::vertex_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from the local vertex
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t vertex,
                       int cell_orientation) const final override
  {
%(tabulate_tensor)s
  }

};
"""

vertex_integral_header = """\
/// This class defines the interface for the tabulation of
/// an expression evaluated at exactly one vertex.

class %(classname)s: public ufc::vertex_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  ~%(classname)s() override;

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override;

  /// Tabulate the tensor for the contribution from the local vertex
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t vertex,
                       int cell_orientation) const final override;

};
"""

vertex_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::vertex_integral()%(initializer_list)s
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
void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
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
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  /// Return the number of cells involved in evaluation of the integral
  std::size_t num_cells() const final override
  {
%(num_cells)s
  }

  /// Tabulate the tensor for the contribution from a custom domain
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       const double * facet_normals,
                       int cell_orientation) const final override
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
  ~%(classname)s() override;

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override;

  /// Return the number of cells involved in evaluation of the integral
  std::size_t num_cells() const final override;

  /// Tabulate the tensor for the contribution from a custom domain
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       const double * facet_normals,
                       int cell_orientation) const final override;

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

/// Tabulate the tensor for the contribution from a custom domain
void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t num_quadrature_points,
                                    const double * quadrature_points,
                                    const double * quadrature_weights,
                                    const double * facet_normals,
                                    int cell_orientation) const
{
%(tabulate_tensor)s
}
"""

cutcell_integral_combined = """\
/// This class defines the interface for the tabulation of the
/// tensor corresponding to the local contribution to a form from
/// the integral over a cut cell defined in terms of a set of
/// quadrature points and weights.

class %(classname)s: public ufc::cutcell_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::cutcell_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from a cutcell domain
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       int cell_orientation) const final override
  {
%(tabulate_tensor)s
  }

};
"""

cutcell_integral_header = """\
/// This class defines the interface for the tabulation of the
/// tensor corresponding to the local contribution to a form from
/// the integral over a cut cell defined in terms of a set of
/// quadrature points and weights.

class %(classname)s: public ufc::cutcell_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  ~%(classname)s() override;

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override;

  /// Tabulate the tensor for the contribution from a cutcell domain
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       int cell_orientation) const final override;

};
"""

cutcell_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::cutcell_integral()%(initializer_list)s
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

/// Tabulate the tensor for the contribution from a cutcell domain
void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t num_quadrature_points,
                                    const double * quadrature_points,
                                    const double * quadrature_weights,
                                    int cell_orientation) const
{
%(tabulate_tensor)s
}
"""

interface_integral_combined = """\
/// This class defines the interface for the tabulation of the
/// tensor corresponding to the local contribution to a form from
/// the integral over an interface cutting through a pair of cells.

class %(classname)s: public ufc::interface_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::interface_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from an interface domain
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       const double * facet_normals,
                       int cell_orientation) const final override
  {
%(tabulate_tensor)s
  }

};
"""

interface_integral_header = """\
/// This class defines the interface for the tabulation of the
/// tensor corresponding to the local contribution to a form from
/// the integral over an interface cutting through a pair of cells.

class %(classname)s: public ufc::interface_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  ~%(classname)s() override;

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override;

  /// Tabulate the tensor for the contribution from an interface domain
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       const double * facet_normals,
                       int cell_orientation) const final override;

};
"""

interface_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::interface_integral()%(initializer_list)s
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

/// Tabulate the tensor for the contribution from an interface domain
void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t num_quadrature_points,
                                    const double * quadrature_points,
                                    const double * quadrature_weights,
                                    const double * facet_normals,
                                    int cell_orientation) const
{
%(tabulate_tensor)s
}
"""

overlap_integral_combined = """\
/// This class defines the interface for the tabulation of the
/// tensor corresponding to the local contribution to a form from
/// the integral over the overlapped portion of a cell defined in
/// terms of a set of quadrature points and weights.

class %(classname)s: public ufc::overlap_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::overlap_integral()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  /// Tabulate the tensor for the contribution from an overlap domain
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       int cell_orientation) const final override
  {
%(tabulate_tensor)s
  }

};
"""

overlap_integral_header = """\
/// This class defines the interface for the tabulation of the
/// tensor corresponding to the local contribution to a form from
/// the integral over the overlapped portion of a cell defined in
/// terms of a set of quadrature points and weights.

class %(classname)s: public ufc::overlap_integral
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  ~%(classname)s() override;

  /// Tabulate which form coefficients are used by this integral
  const std::vector<bool> & enabled_coefficients() const final override;

  /// Tabulate the tensor for the contribution from an overlap domain
  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       int cell_orientation) const final override;

};
"""

overlap_integral_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::overlap_integral()%(initializer_list)s
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

/// Tabulate the tensor for the contribution from an overlap domain
void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t num_quadrature_points,
                                    const double * quadrature_points,
                                    const double * quadrature_weights,
                                    int cell_orientation) const
{
%(tabulate_tensor)s
}
"""
