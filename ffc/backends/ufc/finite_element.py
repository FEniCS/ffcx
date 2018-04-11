# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2017.

finite_element_declaration = """
"""

ufc_finite_element_combined = """
struct {factory_name}
{{
  ufc_finite_element* element = (*ufc_finite_element) malloc(sizeof(*element));
  element->signature = {signature};
  element->ufc_shape = {cell_shape};
  element->topological_dimension = {topological_dimension};
  element->geometric_dimension = {geometric_dimension};
  element->space_dimension = {space_dimension};
  element->value_rank = {value_rank};
  element->value_dimension = {value_dimension}};
  element->value_size = {value_size};
  element->reference_value_rank = {reference_dimension};
  element->reference_value_dimension = {reference_value_dimension};
  element->reference_value_size = {reference_value_size};
  element->degree = {degree};
  element->family = {family};
  element->evaluate_reference_basis = {evaluate_reference_basis};
  element->evaluate_reference_basis_derivatives = {evaluate_reference_basis_derivatives};
  element->transform_reference_basis_derivatives = {transform_reference_basis_derivatives};
  element->map_dofs = {map_dofs};
  element->tabulate_reference_dof_coordinates = {tabulate_reference_dof_coordinates};
  element->num_sub_elements = {num_sub_elements};
  element->create_sub_element = create_sub_element;
  element->create)() = {create};
}};
"""




finite_element_combined = """
class %(classname)s: public ufc::finite_element
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::finite_element()%(initializer_list)s
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

  ufc::shape cell_shape() const final override
  {
%(cell_shape)s
  }

  int64_t topological_dimension() const final override
  {
%(topological_dimension)s
  }

  int64_t geometric_dimension() const final override
  {
%(geometric_dimension)s
  }

  int64_t space_dimension() const final override
  {
%(space_dimension)s
  }

  int64_t value_rank() const final override
  {
%(value_rank)s
  }

  int64_t value_dimension(int64_t i) const final override
  {
%(value_dimension)s
  }

  int64_t value_size() const final override
  {
%(value_size)s
  }

  int64_t reference_value_rank() const final override
  {
%(reference_value_rank)s
  }

  int64_t reference_value_dimension(int64_t i) const final override
  {
%(reference_value_dimension)s
  }

  int64_t reference_value_size() const final override
  {
%(reference_value_size)s
  }

  int64_t degree() const final override
  {
%(degree)s
  }

  const char * family() const final override
  {
%(family)s
  }

  int evaluate_reference_basis(double * reference_values,
                               int64_t num_points,
                               const double * X) const final override
  {
%(evaluate_reference_basis)s
  }

  int evaluate_reference_basis_derivatives(double * reference_values,
                                           int64_t order,
                                           int64_t num_points,
                                           const double * X) const final override
  {
%(evaluate_reference_basis_derivatives)s
  }

  int transform_reference_basis_derivatives(double * values,
                                            int64_t order,
                                            int64_t num_points,
                                            const double * reference_values,
                                            const double * X,
                                            const double * J,
                                            const double * detJ,
                                            const double * K,
                                            int cell_orientation) const final override
  {
%(transform_reference_basis_derivatives)s
  }

  void map_dofs(double * values,
                const double *vals,
                const double * coordinate_dofs,
                int cell_orientation,
                const ufc::coordinate_mapping * cm=NULL
                ) const final override
  {
%(map_dofs)s
  }

  void tabulate_reference_dof_coordinates(double * reference_dof_coordinates) const final override
  {
%(tabulate_reference_dof_coordinates)s
  }

  int64_t num_sub_elements() const final override
  {
%(num_sub_elements)s
  }

  ufc::finite_element * create_sub_element(int64_t i) const final override
  {
%(create_sub_element)s
  }

  ufc::finite_element * create() const final override
  {
%(create)s
  }

};
"""

finite_element_header = """
class %(classname)s: public ufc::finite_element
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const char * signature() const final override;

  ufc::shape cell_shape() const final override;

  int64_t topological_dimension() const final override;

  int64_t geometric_dimension() const final override;

  int64_t space_dimension() const final override;

  int64_t value_rank() const final override;

  int64_t value_dimension(int64_t i) const final override;

  int64_t value_size() const final override;

  int64_t reference_value_rank() const final override;

  int64_t reference_value_dimension(int64_t i) const final override;

  int64_t reference_value_size() const final override;

  int64_t degree() const final override;

  const char * family() const final override;

  int evaluate_reference_basis(double * reference_values,
                               int64_t num_points,
                               const double * X) const final override;

  int evaluate_reference_basis_derivatives(double * reference_values,
                                            int64_t order,
                                            int64_t num_points,
                                            const double * X) const final override;

  int transform_reference_basis_derivatives(double * values,
                                            int64_t order,
                                            int64_t num_points,
                                            const double * reference_values,
                                            const double * X,
                                            const double * J,
                                            const double * detJ,
                                            const double * K,
                                            int cell_orientation) const final override;

  void map_dofs(double * values,
                const double *vals,
                const double * coordinate_dofs,
                int cell_orientation,
                const ufc::coordinate_mapping * cm=NULL
                ) const final override;

  void tabulate_reference_dof_coordinates(double * reference_dof_coordinates) const final override;

  int64_t num_sub_elements() const final override;

  ufc::finite_element * create_sub_element(int64_t i) const final override;

  ufc::finite_element * create() const final override;

};
"""

finite_element_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::finite_element()%(initializer_list)s
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

ufc::shape %(classname)s::cell_shape() const
{
%(cell_shape)s
}

int64_t %(classname)s::topological_dimension() const
{
%(topological_dimension)s
}

int64_t %(classname)s::geometric_dimension() const
{
%(geometric_dimension)s
}

int64_t %(classname)s::space_dimension() const
{
%(space_dimension)s
}

int64_t %(classname)s::value_rank() const
{
%(value_rank)s
}

int64_t %(classname)s::value_dimension(int64_t i) const
{
%(value_dimension)s
}

int64_t %(classname)s::value_size() const
{
%(value_size)s
}

int64_t %(classname)s::reference_value_rank() const
{
%(reference_value_rank)s
}

int64_t %(classname)s::reference_value_dimension(int64_t i) const
{
%(reference_value_dimension)s
}

int64_t %(classname)s::reference_value_size() const
{
%(reference_value_size)s
}

int64_t %(classname)s::degree() const
{
%(degree)s
}

const char * %(classname)s::family() const
{
%(family)s
}

int %(classname)s::evaluate_reference_basis(double * reference_values,
                                            int64_t num_points,
                                            const double * X) const
{
%(evaluate_reference_basis)s
}

int %(classname)s::evaluate_reference_basis_derivatives(double * reference_values,
                                                        int64_t order,
                                                        int64_t num_points,
                                                        const double * X) const
{
%(evaluate_reference_basis_derivatives)s
}

int %(classname)s::transform_reference_basis_derivatives(double * values,
                                                         int64_t order,
                                                         int64_t num_points,
                                                         const double * reference_values,
                                                         const double * X,
                                                         const double * J,
                                                         const double * detJ,
                                                         const double * K,
                                                         int cell_orientation) const
{
%(transform_reference_basis_derivatives)s
}

void %(classname)s::map_dofs(double * values,
                             const double *vals,
                             const double * coordinate_dofs,
                             int cell_orientation,
                             const ufc::coordinate_mapping * cm
                             ) const
{
%(map_dofs)s
}

void %(classname)s::tabulate_reference_dof_coordinates(double * reference_dof_coordinates) const
{
%(tabulate_reference_dof_coordinates)s
}

int64_t %(classname)s::num_sub_elements() const
{
%(num_sub_elements)s
}

ufc::finite_element * %(classname)s::create_sub_element(int64_t i) const
{
%(create_sub_element)s
}

ufc::finite_element * %(classname)s::create() const
{
%(create)s
}
"""
