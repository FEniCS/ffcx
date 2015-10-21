# Code generation format strings for UFC (Unified Form-assembly Code) v. 1.7.0dev.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2015.

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
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Return a string identifying the form
  const char * signature() const final override
  {
%(signature)s
  }


  /// Return the rank of the global tensor (r)
  std::size_t rank() const final override
  {
%(rank)s
  }

  /// Return the number of coefficients (n)
  std::size_t num_coefficients() const final override
  {
%(num_coefficients)s
  }

  /// Return original coefficient position for each coefficient (0 <= i < n)
  std::size_t original_coefficient_position(std::size_t i) const final override
  {
%(original_coefficient_position)s
  }

  /// Create a new finite element for parameterization of coordinates
  ufc::finite_element * create_coordinate_finite_element() const final override
  {
%(create_coordinate_finite_element)s
  }

  /// Create a new dofmap for parameterization of coordinates
  ufc::dofmap * create_coordinate_dofmap() const final override
  {
%(create_coordinate_dofmap)s
   }

  /// Create a new coordinate mapping
  ufc::coordinate_mapping * create_coordinate_mapping() const final override
  {
%(create_coordinate_mapping)s
  }

  /// Create a new finite element for argument function i
  ufc::finite_element * create_finite_element(std::size_t i) const final override
  {
%(create_finite_element)s
  }

  /// Create a new dofmap for argument function i
  ufc::dofmap * create_dofmap(std::size_t i) const final override
  {
%(create_dofmap)s
  }


  /// Return the number of cell domains
  std::size_t max_cell_subdomain_id() const final override
  {
%(max_cell_subdomain_id)s
  }

  /// Return the number of exterior facet domains
  std::size_t max_exterior_facet_subdomain_id() const final override
  {
%(max_exterior_facet_subdomain_id)s
  }

  /// Return the number of interior facet domains
  std::size_t max_interior_facet_subdomain_id() const final override
  {
%(max_interior_facet_subdomain_id)s
  }

  /// Return the number of vertex domains
  std::size_t max_vertex_subdomain_id() const final override
  {
%(max_vertex_subdomain_id)s
  }

  /// Return the number of custom domains
  std::size_t max_custom_subdomain_id() const final override
  {
%(max_custom_subdomain_id)s
  }


  /// Return whether the form has any cell integrals
  bool has_cell_integrals() const final override
  {
%(has_cell_integrals)s
  }

  /// Return whether the form has any exterior facet integrals
  bool has_exterior_facet_integrals() const final override
  {
%(has_exterior_facet_integrals)s
  }

  /// Return whether the form has any interior facet integrals
  bool has_interior_facet_integrals() const final override
  {
%(has_interior_facet_integrals)s
  }

  /// Return whether the form has any vertex integrals
  bool has_vertex_integrals() const final override
  {
%(has_vertex_integrals)s
  }

  /// Return whether the form has any custom integrals
  bool has_custom_integrals() const final override
  {
%(has_custom_integrals)s
  }


  /// Create a new cell integral on sub domain subdomain_id
  ufc::cell_integral * create_cell_integral(std::size_t subdomain_id) const final override
  {
%(create_cell_integral)s
  }

  /// Create a new exterior facet integral on sub domain subdomain_id
  ufc::exterior_facet_integral * create_exterior_facet_integral(std::size_t subdomain_id) const final override
  {
%(create_exterior_facet_integral)s
  }

  /// Create a new interior facet integral on sub domain subdomain_id
  ufc::interior_facet_integral * create_interior_facet_integral(std::size_t subdomain_id) const final override
  {
%(create_interior_facet_integral)s
  }

  /// Create a new vertex integral on sub domain subdomain_id
  ufc::vertex_integral * create_vertex_integral(std::size_t subdomain_id) const final override
  {
%(create_vertex_integral)s
  }

  /// Create a new custom integral on sub domain subdomain_id
  ufc::custom_integral * create_custom_integral(std::size_t subdomain_id) const final override
  {
%(create_custom_integral)s
  }


  /// Create a new cell integral on everywhere else
  ufc::cell_integral * create_default_cell_integral() const final override
  {
%(create_default_cell_integral)s
  }

  /// Create a new exterior facet integral on everywhere else
  ufc::exterior_facet_integral * create_default_exterior_facet_integral() const final override
  {
%(create_default_exterior_facet_integral)s
  }

  /// Create a new interior facet integral on everywhere else
  ufc::interior_facet_integral * create_default_interior_facet_integral() const final override
  {
%(create_default_interior_facet_integral)s
  }

  /// Create a new vertex integral on everywhere else
  ufc::vertex_integral * create_default_vertex_integral() const final override
  {
%(create_default_vertex_integral)s
  }

  /// Create a new custom integral on everywhere else
  ufc::custom_integral * create_default_custom_integral() const final override
  {
%(create_default_custom_integral)s
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
  ~%(classname)s() override;

  /// Return a string identifying the form
  const char * signature() const final override;


  /// Return the rank of the global tensor (r)
  std::size_t rank() const final override;

  /// Return the number of coefficients (n)
  std::size_t num_coefficients() const final override;

  /// Return original coefficient position for each coefficient (0 <= i < n)
  std::size_t original_coefficient_position(std::size_t i) const final override;


  /// Create a new finite element for parameterization of coordinates
  ufc::finite_element * create_coordinate_finite_element() const final override;

  /// Create a new dofmap for parameterization of coordinates
  ufc::dofmap * create_coordinate_dofmap() const final override;

  /// Create a new coordinate mapping
  ufc::coordinate_mapping * create_coordinate_mapping() const final override;

  /// Create a new finite element for argument function i
  ufc::finite_element * create_finite_element(std::size_t i) const final override;

  /// Create a new dofmap for argument function i
  ufc::dofmap * create_dofmap(std::size_t i) const final override;


  /// Return the number of cell domains
  std::size_t max_cell_subdomain_id() const final override;

  /// Return the number of exterior facet domains
  std::size_t max_exterior_facet_subdomain_id() const final override;

  /// Return the number of interior facet domains
  std::size_t max_interior_facet_subdomain_id() const final override;

  /// Return the number of vertex domains
  std::size_t max_vertex_subdomain_id() const final override;

  /// Return the number of custom domains
  std::size_t max_custom_subdomain_id() const final override;


  /// Return whether the form has any cell integrals
  bool has_cell_integrals() const final override;

  /// Return whether the form has any exterior facet integrals
  bool has_exterior_facet_integrals() const final override;

  /// Return whether the form has any interior facet integrals
  bool has_interior_facet_integrals() const final override;

  /// Return whether the form has any vertex integrals
  bool has_vertex_integrals() const final override;

  /// Return whether the form has any custom integrals
  bool has_custom_integrals() const final override;


  /// Create a new cell integral on sub domain i
  ufc::cell_integral * create_cell_integral(std::size_t i) const final override;

  /// Create a new exterior facet integral on sub domain i
  ufc::exterior_facet_integral * create_exterior_facet_integral(std::size_t i) const final override;

  /// Create a new interior facet integral on sub domain i
  ufc::interior_facet_integral * create_interior_facet_integral(std::size_t i) const final override;

  /// Create a new vertex integral on sub domain i
  ufc::vertex_integral * create_vertex_integral(std::size_t i) const final override;

  /// Create a new custom integral on sub domain i
  ufc::custom_integral * create_custom_integral(std::size_t i) const final override;


  /// Create a new cell integral on everywhere else
  ufc::cell_integral * create_default_cell_integral() const final override;

  /// Create a new exterior facet integral on everywhere else
  ufc::exterior_facet_integral * create_default_exterior_facet_integral() const final override;

  /// Create a new interior facet integral on everywhere else
  ufc::interior_facet_integral * create_default_interior_facet_integral() const final override;

  /// Create a new vertex integral on everywhere else
  ufc::vertex_integral * create_default_vertex_integral() const final override;

  /// Create a new custom integral on everywhere else
  ufc::custom_integral * create_default_custom_integral() const final override;
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
const char * %(classname)s::signature() const
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

/// Return original coefficient position for each coefficient (0 <= i < n)
std::size_t %(classname)s::original_coefficient_position(std::size_t i) const
{
%(original_coefficient_position)s
}


/// Create a new finite element for parameterization of coordinates
ufc::finite_element * %(classname)s::create_coordinate_finite_element() const
{
%(create_coordinate_finite_element)s
}

/// Create a new dofmap for parameterization of coordinates
ufc::dofmap * %(classname)s::create_coordinate_dofmap() const
{
%(create_coordinate_dofmap)s
}

/// Create a new coordinate mapping
ufc::coordinate_mapping * %(classname)s::create_coordinate_mapping() const
{
%(create_coordinate_mapping)s
}

/// Create a new finite element for argument function i
ufc::finite_element * %(classname)s::create_finite_element(std::size_t i) const
{
%(create_finite_element)s
}

/// Create a new dofmap for argument function i
ufc::dofmap * %(classname)s::create_dofmap(std::size_t i) const
{
%(create_dofmap)s
}


/// Return the number of cell domains
std::size_t %(classname)s::max_cell_subdomain_id() const
{
%(max_cell_subdomain_id)s
}

/// Return the number of exterior facet domains
std::size_t %(classname)s::max_exterior_facet_subdomain_id() const
{
%(max_exterior_facet_subdomain_id)s
}

/// Return the number of interior facet domains
std::size_t %(classname)s::max_interior_facet_subdomain_id() const
{
%(max_interior_facet_subdomain_id)s
}

/// Return the number of vertex domains
std::size_t %(classname)s::max_vertex_subdomain_id() const
{
%(max_vertex_subdomain_id)s
}

/// Return the number of custom domains
std::size_t %(classname)s::max_custom_subdomain_id() const
{
%(max_custom_subdomain_id)s
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

/// Return whether the form has any vertex integrals
bool %(classname)s::has_vertex_integrals() const
{
%(has_vertex_integrals)s
}

/// Return whether the form has any custom integrals
bool %(classname)s::has_custom_integrals() const
{
%(has_custom_integrals)s
}


/// Create a new cell integral on sub domain subdomain_id
ufc::cell_integral * %(classname)s::create_cell_integral(std::size_t subdomain_id) const
{
%(create_cell_integral)s
}

/// Create a new exterior facet integral on sub domain subdomain_id
ufc::exterior_facet_integral * %(classname)s::create_exterior_facet_integral(std::size_t subdomain_id) const
{
%(create_exterior_facet_integral)s
}

/// Create a new interior facet integral on sub domain subdomain_id
ufc::interior_facet_integral * %(classname)s::create_interior_facet_integral(std::size_t subdomain_id) const
{
%(create_interior_facet_integral)s
}

/// Create a new vertex integral on sub domain subdomain_id
ufc::vertex_integral * %(classname)s::create_vertex_integral(std::size_t subdomain_id) const
{
%(create_vertex_integral)s
}

/// Create a new custom integral on sub domain subdomain_id
ufc::custom_integral * %(classname)s::create_custom_integral(std::size_t subdomain_id) const
{
%(create_custom_integral)s
}


/// Create a new cell integral on everywhere else
ufc::cell_integral * %(classname)s::create_default_cell_integral() const
{
%(create_default_cell_integral)s
}

/// Create a new exterior facet integral on everywhere else
ufc::exterior_facet_integral * %(classname)s::create_default_exterior_facet_integral() const
{
%(create_default_exterior_facet_integral)s
}

/// Create a new interior facet integral on everywhere else
ufc::interior_facet_integral * %(classname)s::create_default_interior_facet_integral() const
{
%(create_default_interior_facet_integral)s
}

/// Create a new vertex integral on everywhere else
ufc::vertex_integral * %(classname)s::create_default_vertex_integral() const
{
%(create_default_vertex_integral)s
}

/// Create a new custom integral on everywhere else
ufc::custom_integral * %(classname)s::create_default_custom_integral() const
{
%(create_default_custom_integral)s
}

"""
