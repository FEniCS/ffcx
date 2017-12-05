# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2017.

form_combined = """
class %(classname)s: public ufc::form
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::form()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const char * signature() const final override
  {
%(signature)s
  }

  std::size_t rank() const final override
  {
%(rank)s
  }

  std::size_t num_coefficients() const final override
  {
%(num_coefficients)s
  }

  std::size_t original_coefficient_position(std::size_t i) const final override
  {
%(original_coefficient_position)s
  }

  ufc::finite_element * create_coordinate_finite_element() const final override
  {
%(create_coordinate_finite_element)s
  }

  ufc::dofmap * create_coordinate_dofmap() const final override
  {
%(create_coordinate_dofmap)s
   }

  ufc::coordinate_mapping * create_coordinate_mapping() const final override
  {
%(create_coordinate_mapping)s
  }

  ufc::finite_element * create_finite_element(std::size_t i) const final override
  {
%(create_finite_element)s
  }

  ufc::dofmap * create_dofmap(std::size_t i) const final override
  {
%(create_dofmap)s
  }

  std::size_t max_cell_subdomain_id() const final override
  {
%(max_cell_subdomain_id)s
  }

  std::size_t max_exterior_facet_subdomain_id() const final override
  {
%(max_exterior_facet_subdomain_id)s
  }

  std::size_t max_interior_facet_subdomain_id() const final override
  {
%(max_interior_facet_subdomain_id)s
  }

  std::size_t max_vertex_subdomain_id() const final override
  {
%(max_vertex_subdomain_id)s
  }

  std::size_t max_custom_subdomain_id() const final override
  {
%(max_custom_subdomain_id)s
  }

  std::size_t max_cutcell_subdomain_id() const final override
  {
%(max_cutcell_subdomain_id)s
  }

  std::size_t max_interface_subdomain_id() const final override
  {
%(max_interface_subdomain_id)s
  }

  std::size_t max_overlap_subdomain_id() const final override
  {
%(max_overlap_subdomain_id)s
  }

  bool has_cell_integrals() const final override
  {
%(has_cell_integrals)s
  }

  bool has_exterior_facet_integrals() const final override
  {
%(has_exterior_facet_integrals)s
  }

  bool has_interior_facet_integrals() const final override
  {
%(has_interior_facet_integrals)s
  }

  bool has_vertex_integrals() const final override
  {
%(has_vertex_integrals)s
  }

  bool has_custom_integrals() const final override
  {
%(has_custom_integrals)s
  }

  bool has_cutcell_integrals() const final override
  {
%(has_cutcell_integrals)s
  }

  bool has_interface_integrals() const final override
  {
%(has_interface_integrals)s
  }

  bool has_overlap_integrals() const final override
  {
%(has_overlap_integrals)s
  }

  ufc::cell_integral * create_cell_integral(std::size_t subdomain_id) const final override
  {
%(create_cell_integral)s
  }

  ufc::exterior_facet_integral * create_exterior_facet_integral(std::size_t subdomain_id) const final override
  {
%(create_exterior_facet_integral)s
  }

  ufc::interior_facet_integral * create_interior_facet_integral(std::size_t subdomain_id) const final override
  {
%(create_interior_facet_integral)s
  }

  ufc::vertex_integral * create_vertex_integral(std::size_t subdomain_id) const final override
  {
%(create_vertex_integral)s
  }

  ufc::custom_integral * create_custom_integral(std::size_t subdomain_id) const final override
  {
%(create_custom_integral)s
  }

  ufc::cutcell_integral * create_cutcell_integral(std::size_t subdomain_id) const final override
  {
%(create_cutcell_integral)s
  }

  ufc::interface_integral * create_interface_integral(std::size_t subdomain_id) const final override
  {
%(create_interface_integral)s
  }

  ufc::overlap_integral * create_overlap_integral(std::size_t subdomain_id) const final override
  {
%(create_overlap_integral)s
  }

  ufc::cell_integral * create_default_cell_integral() const final override
  {
%(create_default_cell_integral)s
  }

  ufc::exterior_facet_integral * create_default_exterior_facet_integral() const final override
  {
%(create_default_exterior_facet_integral)s
  }

  ufc::interior_facet_integral * create_default_interior_facet_integral() const final override
  {
%(create_default_interior_facet_integral)s
  }

  ufc::vertex_integral * create_default_vertex_integral() const final override
  {
%(create_default_vertex_integral)s
  }

  ufc::custom_integral * create_default_custom_integral() const final override
  {
%(create_default_custom_integral)s
  }

  ufc::cutcell_integral * create_default_cutcell_integral() const final override
  {
%(create_default_cutcell_integral)s
  }

  ufc::interface_integral * create_default_interface_integral() const final override
  {
%(create_default_interface_integral)s
  }

  ufc::overlap_integral * create_default_overlap_integral() const final override
  {
%(create_default_overlap_integral)s
  }

};
"""

form_header = """
class %(classname)s: public ufc::form
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const char * signature() const final override;

  std::size_t rank() const final override;

  std::size_t num_coefficients() const final override;

  std::size_t original_coefficient_position(std::size_t i) const final override;

  ufc::finite_element * create_coordinate_finite_element() const final override;

  ufc::dofmap * create_coordinate_dofmap() const final override;

  ufc::coordinate_mapping * create_coordinate_mapping() const final override;

  ufc::finite_element * create_finite_element(std::size_t i) const final override;

  ufc::dofmap * create_dofmap(std::size_t i) const final override;

  std::size_t max_cell_subdomain_id() const final override;

  std::size_t max_exterior_facet_subdomain_id() const final override;

  std::size_t max_interior_facet_subdomain_id() const final override;

  std::size_t max_vertex_subdomain_id() const final override;

  std::size_t max_custom_subdomain_id() const final override;

  std::size_t max_cutcell_subdomain_id() const final override;

  std::size_t max_interface_subdomain_id() const final override;

  std::size_t max_overlap_subdomain_id() const final override;

  bool has_cell_integrals() const final override;

  bool has_exterior_facet_integrals() const final override;

  bool has_interior_facet_integrals() const final override;

  bool has_vertex_integrals() const final override;

  bool has_custom_integrals() const final override;

  bool has_cutcell_integrals() const final override;

  bool has_interface_integrals() const final override;

  bool has_overlap_integrals() const final override;

  ufc::cell_integral * create_cell_integral(std::size_t i) const final override;

  ufc::exterior_facet_integral * create_exterior_facet_integral(std::size_t i) const final override;

  ufc::interior_facet_integral * create_interior_facet_integral(std::size_t i) const final override;

  ufc::vertex_integral * create_vertex_integral(std::size_t i) const final override;

  ufc::custom_integral * create_custom_integral(std::size_t i) const final override;

  ufc::cutcell_integral * create_cutcell_integral(std::size_t i) const final override;

  ufc::interface_integral * create_interface_integral(std::size_t i) const final override;

  ufc::overlap_integral * create_overlap_integral(std::size_t i) const final override;

  ufc::cell_integral * create_default_cell_integral() const final override;

  ufc::exterior_facet_integral * create_default_exterior_facet_integral() const final override;

  ufc::interior_facet_integral * create_default_interior_facet_integral() const final override;

  ufc::vertex_integral * create_default_vertex_integral() const final override;

  ufc::custom_integral * create_default_custom_integral() const final override;

  ufc::cutcell_integral * create_default_cutcell_integral() const final override;

  ufc::interface_integral * create_default_interface_integral() const final override;

  ufc::overlap_integral * create_default_overlap_integral() const final override;

};
"""

form_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::form()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const char * %(classname)s::signature() const
{
%(signature)s
}

std::size_t %(classname)s::rank() const
{
%(rank)s
}

std::size_t %(classname)s::num_coefficients() const
{
%(num_coefficients)s
}

std::size_t %(classname)s::original_coefficient_position(std::size_t i) const
{
%(original_coefficient_position)s
}

ufc::finite_element * %(classname)s::create_coordinate_finite_element() const
{
%(create_coordinate_finite_element)s
}

ufc::dofmap * %(classname)s::create_coordinate_dofmap() const
{
%(create_coordinate_dofmap)s
}

ufc::coordinate_mapping * %(classname)s::create_coordinate_mapping() const
{
%(create_coordinate_mapping)s
}

ufc::finite_element * %(classname)s::create_finite_element(std::size_t i) const
{
%(create_finite_element)s
}

ufc::dofmap * %(classname)s::create_dofmap(std::size_t i) const
{
%(create_dofmap)s
}

std::size_t %(classname)s::max_cell_subdomain_id() const
{
%(max_cell_subdomain_id)s
}

std::size_t %(classname)s::max_exterior_facet_subdomain_id() const
{
%(max_exterior_facet_subdomain_id)s
}

std::size_t %(classname)s::max_interior_facet_subdomain_id() const
{
%(max_interior_facet_subdomain_id)s
}

std::size_t %(classname)s::max_vertex_subdomain_id() const
{
%(max_vertex_subdomain_id)s
}

std::size_t %(classname)s::max_custom_subdomain_id() const
{
%(max_custom_subdomain_id)s
}

std::size_t %(classname)s::max_cutcell_subdomain_id() const
{
%(max_cutcell_subdomain_id)s
}

std::size_t %(classname)s::max_interface_subdomain_id() const
{
%(max_interface_subdomain_id)s
}

std::size_t %(classname)s::max_overlap_subdomain_id() const
{
%(max_overlap_subdomain_id)s
}

bool %(classname)s::has_cell_integrals() const
{
%(has_cell_integrals)s
}

bool %(classname)s::has_exterior_facet_integrals() const
{
%(has_exterior_facet_integrals)s
}

bool %(classname)s::has_interior_facet_integrals() const
{
%(has_interior_facet_integrals)s
}

bool %(classname)s::has_vertex_integrals() const
{
%(has_vertex_integrals)s
}

bool %(classname)s::has_custom_integrals() const
{
%(has_custom_integrals)s
}

bool %(classname)s::has_cutcell_integrals() const
{
%(has_cutcell_integrals)s
}

bool %(classname)s::has_interface_integrals() const
{
%(has_interface_integrals)s
}

bool %(classname)s::has_overlap_integrals() const
{
%(has_overlap_integrals)s
}

ufc::cell_integral * %(classname)s::create_cell_integral(std::size_t subdomain_id) const
{
%(create_cell_integral)s
}

ufc::exterior_facet_integral * %(classname)s::create_exterior_facet_integral(std::size_t subdomain_id) const
{
%(create_exterior_facet_integral)s
}

ufc::interior_facet_integral * %(classname)s::create_interior_facet_integral(std::size_t subdomain_id) const
{
%(create_interior_facet_integral)s
}

ufc::vertex_integral * %(classname)s::create_vertex_integral(std::size_t subdomain_id) const
{
%(create_vertex_integral)s
}

ufc::custom_integral * %(classname)s::create_custom_integral(std::size_t subdomain_id) const
{
%(create_custom_integral)s
}

ufc::cutcell_integral * %(classname)s::create_cutcell_integral(std::size_t subdomain_id) const
{
%(create_cutcell_integral)s
}

ufc::interface_integral * %(classname)s::create_interface_integral(std::size_t subdomain_id) const
{
%(create_interface_integral)s
}

ufc::overlap_integral * %(classname)s::create_overlap_integral(std::size_t subdomain_id) const
{
%(create_overlap_integral)s
}

ufc::cell_integral * %(classname)s::create_default_cell_integral() const
{
%(create_default_cell_integral)s
}

ufc::exterior_facet_integral * %(classname)s::create_default_exterior_facet_integral() const
{
%(create_default_exterior_facet_integral)s
}

ufc::interior_facet_integral * %(classname)s::create_default_interior_facet_integral() const
{
%(create_default_interior_facet_integral)s
}

ufc::vertex_integral * %(classname)s::create_default_vertex_integral() const
{
%(create_default_vertex_integral)s
}

ufc::custom_integral * %(classname)s::create_default_custom_integral() const
{
%(create_default_custom_integral)s
}

ufc::cutcell_integral * %(classname)s::create_default_cutcell_integral() const
{
%(create_default_cutcell_integral)s
}

ufc::interface_integral * %(classname)s::create_default_interface_integral() const
{
%(create_default_interface_integral)s
}

ufc::overlap_integral * %(classname)s::create_default_overlap_integral() const
{
%(create_default_overlap_integral)s
}
"""
