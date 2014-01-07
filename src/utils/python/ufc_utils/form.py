# Code generation format strings for UFC (Unified Form-assembly Code) v. 2.3.0.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2014.

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
  virtual std::size_t rank() const
  {
%(rank)s
  }

  /// Return the number of coefficients (n)
  virtual std::size_t num_coefficients() const
  {
%(num_coefficients)s
  }

  /// Return the number of cell domains
  virtual std::size_t num_cell_domains() const
  {
%(num_cell_domains)s
  }

  /// Return the number of exterior facet domains
  virtual std::size_t num_exterior_facet_domains() const
  {
%(num_exterior_facet_domains)s
  }

  /// Return the number of interior facet domains
  virtual std::size_t num_interior_facet_domains() const
  {
%(num_interior_facet_domains)s
  }

  /// Return the number of point domains
  virtual std::size_t num_point_domains() const
  {
%(num_point_domains)s
  }

  /// Return whether the form has any cell integrals
  virtual bool has_cell_integrals() const
  {
%(has_cell_integrals)s
  }

  /// Return whether the form has any exterior facet integrals
  virtual bool has_exterior_facet_integrals() const
  {
%(has_exterior_facet_integrals)s
  }

  /// Return whether the form has any interior facet integrals
  virtual bool has_interior_facet_integrals() const
  {
%(has_interior_facet_integrals)s
  }

  /// Return whether the form has any point integrals
  virtual bool has_point_integrals() const
  {
%(has_point_integrals)s
  }

  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(std::size_t i) const
  {
%(create_finite_element)s
  }

  /// Create a new dofmap for argument function i
  virtual ufc::dofmap* create_dofmap(std::size_t i) const
  {
%(create_dofmap)s
  }

  /// Create a new cell integral on sub domain i
  virtual ufc::cell_integral* create_cell_integral(std::size_t i) const
  {
%(create_cell_integral)s
  }

  /// Create a new exterior facet integral on sub domain i
  virtual ufc::exterior_facet_integral* create_exterior_facet_integral(std::size_t i) const
  {
%(create_exterior_facet_integral)s
  }

  /// Create a new interior facet integral on sub domain i
  virtual ufc::interior_facet_integral* create_interior_facet_integral(std::size_t i) const
  {
%(create_interior_facet_integral)s
  }

  /// Create a new point integral on sub domain i
  virtual ufc::point_integral* create_point_integral(std::size_t i) const
  {
%(create_point_integral)s
  }

  /// Create a new cell integral on everywhere else
  virtual ufc::cell_integral* create_default_cell_integral() const
  {
%(create_default_cell_integral)s
  }

  /// Create a new exterior facet integral on everywhere else
  virtual ufc::exterior_facet_integral* create_default_exterior_facet_integral() const
  {
%(create_default_exterior_facet_integral)s
  }

  /// Create a new interior facet integral on everywhere else
  virtual ufc::interior_facet_integral* create_default_interior_facet_integral() const
  {
%(create_default_interior_facet_integral)s
  }

  /// Create a new point integral on everywhere else
  virtual ufc::point_integral* create_default_point_integral() const
  {
%(create_default_point_integral)s
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
  virtual std::size_t rank() const;

  /// Return the number of coefficients (n)
  virtual std::size_t num_coefficients() const;

  /// Return the number of cell domains
  virtual std::size_t num_cell_domains() const;

  /// Return the number of exterior facet domains
  virtual std::size_t num_exterior_facet_domains() const;

  /// Return the number of interior facet domains
  virtual std::size_t num_interior_facet_domains() const;

  /// Return the number of point domains
  virtual std::size_t num_point_domains() const;

  /// Return whether the form has any cell integrals
  virtual bool has_cell_integrals() const;

  /// Return whether the form has any exterior facet integrals
  virtual bool has_exterior_facet_integrals() const;

  /// Return whether the form has any interior facet integrals
  virtual bool has_interior_facet_integrals() const;

  /// Return whether the form has any point integrals
  virtual bool has_point_integrals() const;

  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(std::size_t i) const;

  /// Create a new dofmap for argument function i
  virtual ufc::dofmap* create_dofmap(std::size_t i) const;

  /// Create a new cell integral on sub domain i
  virtual ufc::cell_integral* create_cell_integral(std::size_t i) const;

  /// Create a new exterior facet integral on sub domain i
  virtual ufc::exterior_facet_integral* create_exterior_facet_integral(std::size_t i) const;

  /// Create a new interior facet integral on sub domain i
  virtual ufc::interior_facet_integral* create_interior_facet_integral(std::size_t i) const;

  /// Create a new point integral on sub domain i
  virtual ufc::point_integral* create_point_integral(std::size_t i) const;

  /// Create a new cell integral on everywhere else
  virtual ufc::cell_integral* create_default_cell_integral() const;

  /// Create a new exterior facet integral on everywhere else
  virtual ufc::exterior_facet_integral* create_default_exterior_facet_integral() const;

  /// Create a new interior facet integral on everywhere else
  virtual ufc::interior_facet_integral* create_default_interior_facet_integral() const;

  /// Create a new point integral on everywhere else
  virtual ufc::point_integral* create_default_point_integral() const;
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
std::size_t %(classname)s::rank() const
{
%(rank)s
}

/// Return the number of coefficients (n)
std::size_t %(classname)s::num_coefficients() const
{
%(num_coefficients)s
}

/// Return the number of cell domains
std::size_t %(classname)s::num_cell_domains() const
{
%(num_cell_domains)s
}

/// Return the number of exterior facet domains
std::size_t %(classname)s::num_exterior_facet_domains() const
{
%(num_exterior_facet_domains)s
}

/// Return the number of interior facet domains
std::size_t %(classname)s::num_interior_facet_domains() const
{
%(num_interior_facet_domains)s
}

/// Return the number of point domains
std::size_t %(classname)s::num_point_domains() const
{
%(num_point_domains)s
}

/// Return whether the form has any cell integrals
bool %(classname)s::has_cell_integrals() const
{
%(has_cell_integrals)s
}

/// Return whether the form has any exterior facet integrals
bool %(classname)s::has_exterior_facet_integrals() const
{
%(has_exterior_facet_integrals)s
}

/// Return whether the form has any interior facet integrals
bool %(classname)s::has_interior_facet_integrals() const
{
%(has_interior_facet_integrals)s
}

/// Return whether the form has any point integrals
bool %(classname)s::has_point_integrals() const
{
%(has_point_integrals)s
}

/// Create a new finite element for argument function i
ufc::finite_element* %(classname)s::create_finite_element(std::size_t i) const
{
%(create_finite_element)s
}

/// Create a new dofmap for argument function i
ufc::dofmap* %(classname)s::create_dofmap(std::size_t i) const
{
%(create_dofmap)s
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* %(classname)s::create_cell_integral(std::size_t i) const
{
%(create_cell_integral)s
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* %(classname)s::create_exterior_facet_integral(std::size_t i) const
{
%(create_exterior_facet_integral)s
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* %(classname)s::create_interior_facet_integral(std::size_t i) const
{
%(create_interior_facet_integral)s
}

/// Create a new point integral on sub domain i
ufc::point_integral* %(classname)s::create_point_integral(std::size_t i) const
{
%(create_point_integral)s
}

/// Create a new cell integral on everywhere else
ufc::cell_integral* %(classname)s::create_default_cell_integral() const
{
%(create_default_cell_integral)s
}

/// Create a new exterior facet integral on everywhere else
ufc::exterior_facet_integral* %(classname)s::create_default_exterior_facet_integral() const
{
%(create_default_exterior_facet_integral)s
}

/// Create a new interior facet integral on everywhere else
ufc::interior_facet_integral* %(classname)s::create_default_interior_facet_integral() const
{
%(create_default_interior_facet_integral)s
}

/// Create a new point integral on everywhere else
ufc::point_integral* %(classname)s::create_default_point_integral() const
{
%(create_default_point_integral)s
}

"""
