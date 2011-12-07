# Code generation format strings for UFC (Unified Form-assembly Code) v. 2.0.5.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2011.

form_combined = """\
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

class %(classname)s: public ufc::form
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::form()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Return a string identifying the form
  virtual const char* signature() const
  {
%(signature)s
  }

  /// Return the rank of the global tensor (r)
  virtual unsigned int rank() const
  {
%(rank)s
  }

  /// Return the number of coefficients (n)
  virtual unsigned int num_coefficients() const
  {
%(num_coefficients)s
  }

  /// Return the number of cell domains
  virtual unsigned int num_cell_domains() const
  {
%(num_cell_domains)s
  }

  /// Return the number of exterior facet domains
  virtual unsigned int num_exterior_facet_domains() const
  {
%(num_exterior_facet_domains)s
  }

  /// Return the number of interior facet domains
  virtual unsigned int num_interior_facet_domains() const
  {
%(num_interior_facet_domains)s
  }

  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(unsigned int i) const
  {
%(create_finite_element)s
  }

  /// Create a new dofmap for argument function i
  virtual ufc::dofmap* create_dofmap(unsigned int i) const
  {
%(create_dofmap)s
  }

  /// Create a new cell integral on sub domain i
  virtual ufc::cell_integral* create_cell_integral(unsigned int i) const
  {
%(create_cell_integral)s
  }

  /// Create a new exterior facet integral on sub domain i
  virtual ufc::exterior_facet_integral* create_exterior_facet_integral(unsigned int i) const
  {
%(create_exterior_facet_integral)s
  }

  /// Create a new interior facet integral on sub domain i
  virtual ufc::interior_facet_integral* create_interior_facet_integral(unsigned int i) const
  {
%(create_interior_facet_integral)s
  }

};
"""

form_header = """\
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

class %(classname)s: public ufc::form
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  virtual ~%(classname)s();

  /// Return a string identifying the form
  virtual const char* signature() const;

  /// Return the rank of the global tensor (r)
  virtual unsigned int rank() const;

  /// Return the number of coefficients (n)
  virtual unsigned int num_coefficients() const;

  /// Return the number of cell domains
  virtual unsigned int num_cell_domains() const;

  /// Return the number of exterior facet domains
  virtual unsigned int num_exterior_facet_domains() const;

  /// Return the number of interior facet domains
  virtual unsigned int num_interior_facet_domains() const;

  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(unsigned int i) const;

  /// Create a new dofmap for argument function i
  virtual ufc::dofmap* create_dofmap(unsigned int i) const;

  /// Create a new cell integral on sub domain i
  virtual ufc::cell_integral* create_cell_integral(unsigned int i) const;

  /// Create a new exterior facet integral on sub domain i
  virtual ufc::exterior_facet_integral* create_exterior_facet_integral(unsigned int i) const;

  /// Create a new interior facet integral on sub domain i
  virtual ufc::interior_facet_integral* create_interior_facet_integral(unsigned int i) const;

};
"""

form_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::form()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Return a string identifying the form
const char* %(classname)s::signature() const
{
%(signature)s
}

/// Return the rank of the global tensor (r)
unsigned int %(classname)s::rank() const
{
%(rank)s
}

/// Return the number of coefficients (n)
unsigned int %(classname)s::num_coefficients() const
{
%(num_coefficients)s
}

/// Return the number of cell domains
unsigned int %(classname)s::num_cell_domains() const
{
%(num_cell_domains)s
}

/// Return the number of exterior facet domains
unsigned int %(classname)s::num_exterior_facet_domains() const
{
%(num_exterior_facet_domains)s
}

/// Return the number of interior facet domains
unsigned int %(classname)s::num_interior_facet_domains() const
{
%(num_interior_facet_domains)s
}

/// Create a new finite element for argument function i
ufc::finite_element* %(classname)s::create_finite_element(unsigned int i) const
{
%(create_finite_element)s
}

/// Create a new dofmap for argument function i
ufc::dofmap* %(classname)s::create_dofmap(unsigned int i) const
{
%(create_dofmap)s
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* %(classname)s::create_cell_integral(unsigned int i) const
{
%(create_cell_integral)s
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* %(classname)s::create_exterior_facet_integral(unsigned int i) const
{
%(create_exterior_facet_integral)s
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* %(classname)s::create_interior_facet_integral(unsigned int i) const
{
%(create_interior_facet_integral)s
}
"""
