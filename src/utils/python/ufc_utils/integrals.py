# Code generation format strings for UFC (Unified Form-assembly Code) v. 2.0.5.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2011.

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

  /// Tabulate the tensor for the contribution from a local cell
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c) const
  {
%(tabulate_tensor)s
  }

  /// Tabulate the tensor for the contribution from a local cell
  /// using the specified reference cell quadrature points/weights
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int num_quadrature_points,
                               const double * const * quadrature_points,
                               const double* quadrature_weights) const
  {
%(tabulate_tensor_quadrature)s
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

  /// Tabulate the tensor for the contribution from a local cell
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c) const;

  /// Tabulate the tensor for the contribution from a local cell
  /// using the specified reference cell quadrature points/weights
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int num_quadrature_points,
                               const double * const * quadrature_points,
                               const double* quadrature_weights) const;

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

/// Tabulate the tensor for the contribution from a local cell
void %(classname)s::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
%(tabulate_tensor)s
}

/// Tabulate the tensor for the contribution from a local cell
/// using the specified reference cell quadrature points/weights
void %(classname)s::tabulate_tensor(double* A,
                     const double * const * w,
                     const ufc::cell& c,
                     unsigned int num_quadrature_points,
                     const double * const * quadrature_points,
                     const double* quadrature_weights) const
{
%(tabulate_tensor_quadrature)s
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

  /// Tabulate the tensor for the contribution from a local exterior facet
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int facet) const
  {
%(tabulate_tensor)s
  }

  /// Tabulate the tensor for the contribution from a local exterior facet
  /// using the specified reference cell quadrature points/weights
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int num_quadrature_points,
                               const double * const * quadrature_points,
                               const double* quadrature_weights) const
  {
%(tabulate_tensor_quadrature)s
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

  /// Tabulate the tensor for the contribution from a local exterior facet
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int facet) const;

  /// Tabulate the tensor for the contribution from a local exterior facet
  /// using the specified reference cell quadrature points/weights
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int num_quadrature_points,
                               const double * const * quadrature_points,
                               const double* quadrature_weights) const;

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

/// Tabulate the tensor for the contribution from a local exterior facet
void %(classname)s::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
%(tabulate_tensor)s
}

/// Tabulate the tensor for the contribution from a local exterior facet
/// using the specified reference cell quadrature points/weights
void %(classname)s::tabulate_tensor(double* A,
                     const double * const * w,
                     const ufc::cell& c,
                     unsigned int num_quadrature_points,
                     const double * const * quadrature_points,
                     const double* quadrature_weights) const
{
%(tabulate_tensor_quadrature)s
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

  /// Tabulate the tensor for the contribution from a local interior facet
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c0,
                               const ufc::cell& c1,
                               unsigned int facet0,
                               unsigned int facet1) const
  {
%(tabulate_tensor)s
  }

  /// Tabulate the tensor for the contribution from a local interior facet
  /// using the specified reference cell quadrature points/weights
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int num_quadrature_points,
                               const double * const * quadrature_points,
                               const double* quadrature_weights) const
  {
%(tabulate_tensor_quadrature)s
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

  /// Tabulate the tensor for the contribution from a local interior facet
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c0,
                               const ufc::cell& c1,
                               unsigned int facet0,
                               unsigned int facet1) const;

  /// Tabulate the tensor for the contribution from a local interior facet
  /// using the specified reference cell quadrature points/weights
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int num_quadrature_points,
                               const double * const * quadrature_points,
                               const double* quadrature_weights) const;

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

/// Tabulate the tensor for the contribution from a local interior facet
void %(classname)s::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c0,
                                    const ufc::cell& c1,
                                    unsigned int facet0,
                                    unsigned int facet1) const
{
%(tabulate_tensor)s
}

/// Tabulate the tensor for the contribution from a local interior facet
/// using the specified reference cell quadrature points/weights
void %(classname)s::tabulate_tensor(double* A,
                     const double * const * w,
                     const ufc::cell& c,
                     unsigned int num_quadrature_points,
                     const double * const * quadrature_points,
                     const double* quadrature_weights) const
{
%(tabulate_tensor_quadrature)s
}
"""
