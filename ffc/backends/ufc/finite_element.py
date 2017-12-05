# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2017.

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

  std::size_t topological_dimension() const final override
  {
%(topological_dimension)s
  }

  std::size_t geometric_dimension() const final override
  {
%(geometric_dimension)s
  }

  std::size_t space_dimension() const final override
  {
%(space_dimension)s
  }

  std::size_t value_rank() const final override
  {
%(value_rank)s
  }

  std::size_t value_dimension(std::size_t i) const final override
  {
%(value_dimension)s
  }

  std::size_t value_size() const final override
  {
%(value_size)s
  }

  std::size_t reference_value_rank() const final override
  {
%(reference_value_rank)s
  }

  std::size_t reference_value_dimension(std::size_t i) const final override
  {
%(reference_value_dimension)s
  }

  std::size_t reference_value_size() const final override
  {
%(reference_value_size)s
  }

  std::size_t degree() const final override
  {
%(degree)s
  }

  const char * family() const final override
  {
%(family)s
  }

  void evaluate_reference_basis(double * reference_values,
                                std::size_t num_points,
                                const double * X) const final override
  {
%(evaluate_reference_basis)s
  }

  void evaluate_reference_basis_derivatives(double * reference_values,
                                            std::size_t order,
                                            std::size_t num_points,
                                            const double * X) const final override
  {
%(evaluate_reference_basis_derivatives)s
  }

  void transform_reference_basis_derivatives(double * values,
                                             std::size_t order,
                                             std::size_t num_points,
                                             const double * reference_values,
                                             const double * X,
                                             const double * J,
                                             const double * detJ,
                                             const double * K,
                                             int cell_orientation) const final override
  {
%(transform_reference_basis_derivatives)s
  }

  void evaluate_basis(std::size_t i,
                      double * values,
                      const double * x,
                      const double * coordinate_dofs,
                      int cell_orientation,
                      const ufc::coordinate_mapping * cm=nullptr
                      ) const final override
  {
%(evaluate_basis)s
  }

  void evaluate_basis_all(double * values,
                          const double * x,
                          const double * coordinate_dofs,
                          int cell_orientation,
                          const ufc::coordinate_mapping * cm=nullptr
                          ) const final override
  {
%(evaluate_basis_all)s
  }

  void evaluate_basis_derivatives(std::size_t i,
                                  std::size_t n,
                                  double * values,
                                  const double * x,
                                  const double * coordinate_dofs,
                                  int cell_orientation,
                                  const ufc::coordinate_mapping * cm=nullptr
                                  ) const final override
  {
%(evaluate_basis_derivatives)s
  }

  void evaluate_basis_derivatives_all(std::size_t n,
                                      double * values,
                                      const double * x,
                                      const double * coordinate_dofs,
                                      int cell_orientation,
                                      const ufc::coordinate_mapping * cm=nullptr
                                      ) const final override
  {
%(evaluate_basis_derivatives_all)s
  }

  double evaluate_dof(std::size_t i,
                      const ufc::function& f,
                      const double * coordinate_dofs,
                      int cell_orientation,
                      const ufc::cell& c,
                      const ufc::coordinate_mapping * cm=nullptr
                      ) const final override
  {
%(evaluate_dof)s
  }

  void evaluate_dofs(double * values,
                             const ufc::function& f,
                             const double * coordinate_dofs,
                             int cell_orientation,
                             const ufc::cell& c,
                             const ufc::coordinate_mapping * cm=nullptr
                             ) const final override
  {
%(evaluate_dofs)s
  }

  void interpolate_vertex_values(double * vertex_values,
                                 const double * dof_values,
                                 const double * coordinate_dofs,
                                 int cell_orientation,
                                 const ufc::coordinate_mapping * cm=nullptr
                                 ) const final override
  {
%(interpolate_vertex_values)s
  }

  void tabulate_dof_coordinates(double * dof_coordinates,
                                const double * coordinate_dofs,
                                const ufc::coordinate_mapping * cm=nullptr
                                ) const final override
  {
%(tabulate_dof_coordinates)s
  }

  void tabulate_reference_dof_coordinates(double * reference_dof_coordinates) const final override
  {
%(tabulate_reference_dof_coordinates)s
  }

  std::size_t num_sub_elements() const final override
  {
%(num_sub_elements)s
  }

  ufc::finite_element * create_sub_element(std::size_t i) const final override
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

  std::size_t topological_dimension() const final override;

  std::size_t geometric_dimension() const final override;

  std::size_t space_dimension() const final override;

  std::size_t value_rank() const final override;

  std::size_t value_dimension(std::size_t i) const final override;

  std::size_t value_size() const final override;

  std::size_t reference_value_rank() const final override;

  std::size_t reference_value_dimension(std::size_t i) const final override;

  std::size_t reference_value_size() const final override;

  std::size_t degree() const final override;

  const char * family() const final override;

  void evaluate_reference_basis(double * reference_values,
                                std::size_t num_points,
                                const double * X) const final override;

  void evaluate_reference_basis_derivatives(double * reference_values,
                                            std::size_t order,
                                            std::size_t num_points,
                                            const double * X) const final override;

  void transform_reference_basis_derivatives(double * values,
                                             std::size_t order,
                                             std::size_t num_points,
                                             const double * reference_values,
                                             const double * X,
                                             const double * J,
                                             const double * detJ,
                                             const double * K,
                                             int cell_orientation) const final override;

  void evaluate_basis(std::size_t i,
                      double * values,
                      const double * x,
                      const double * coordinate_dofs,
                      int cell_orientation,
                      const ufc::coordinate_mapping * cm=nullptr
                      ) const final override;

  void evaluate_basis_all(double * values,
                          const double * x,
                          const double * coordinate_dofs,
                          int cell_orientation,
                          const ufc::coordinate_mapping * cm=nullptr
                          ) const final override;

  void evaluate_basis_derivatives(std::size_t i,
                                  std::size_t n,
                                  double * values,
                                  const double * x,
                                  const double * coordinate_dofs,
                                  int cell_orientation,
                                  const ufc::coordinate_mapping * cm=nullptr
                                  ) const final override;

  void evaluate_basis_derivatives_all(std::size_t n,
                                      double * values,
                                      const double * x,
                                      const double * coordinate_dofs,
                                      int cell_orientation,
                                      const ufc::coordinate_mapping * cm=nullptr
                                      ) const final override;

  double evaluate_dof(std::size_t i,
                      const ufc::function& f,
                      const double * coordinate_dofs,
                      int cell_orientation,
                      const ufc::cell& c,
                      const ufc::coordinate_mapping * cm=nullptr
                      ) const final override;

  void evaluate_dofs(double * values,
                     const ufc::function& f,
                     const double * coordinate_dofs,
                     int cell_orientation,
                     const ufc::cell& c,
                     const ufc::coordinate_mapping * cm=nullptr
                     ) const final override;

  void interpolate_vertex_values(double * vertex_values,
                                 const double * dof_values,
                                 const double * coordinate_dofs,
                                 int cell_orientation,
                                 const ufc::coordinate_mapping * cm=nullptr
                                 ) const final override;

  void tabulate_dof_coordinates(double * dof_coordinates,
                                const double * coordinate_dofs,
                                const ufc::coordinate_mapping * cm=nullptr
                                ) const final override;

  void tabulate_reference_dof_coordinates(double * reference_dof_coordinates) const final override;

  std::size_t num_sub_elements() const final override;

  ufc::finite_element * create_sub_element(std::size_t i) const final override;

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

std::size_t %(classname)s::topological_dimension() const
{
%(topological_dimension)s
}

std::size_t %(classname)s::geometric_dimension() const
{
%(geometric_dimension)s
}

std::size_t %(classname)s::space_dimension() const
{
%(space_dimension)s
}

std::size_t %(classname)s::value_rank() const
{
%(value_rank)s
}

std::size_t %(classname)s::value_dimension(std::size_t i) const
{
%(value_dimension)s
}

std::size_t %(classname)s::value_size() const
{
%(value_size)s
}

std::size_t %(classname)s::reference_value_rank() const
{
%(reference_value_rank)s
}

std::size_t %(classname)s::reference_value_dimension(std::size_t i) const
{
%(reference_value_dimension)s
}

std::size_t %(classname)s::reference_value_size() const
{
%(reference_value_size)s
}

std::size_t %(classname)s::degree() const
{
%(degree)s
}

const char * %(classname)s::family() const
{
%(family)s
}

void %(classname)s::evaluate_reference_basis(double * reference_values,
                                             std::size_t num_points,
                                             const double * X) const
{
%(evaluate_reference_basis)s
}

void %(classname)s::evaluate_reference_basis_derivatives(double * reference_values,
                                                         std::size_t order,
                                                         std::size_t num_points,
                                                         const double * X) const
{
%(evaluate_reference_basis_derivatives)s
}

void %(classname)s::transform_reference_basis_derivatives(double * values,
                                                          std::size_t order,
                                                          std::size_t num_points,
                                                          const double * reference_values,
                                                          const double * X,
                                                          const double * J,
                                                          const double * detJ,
                                                          const double * K,
                                                          int cell_orientation) const
{
%(transform_reference_basis_derivatives)s
}

void %(classname)s::evaluate_basis(std::size_t i,
                                   double * values,
                                   const double * x,
                                   const double * coordinate_dofs,
                                   int cell_orientation,
                                   const ufc::coordinate_mapping * cm
                                   ) const
{
%(evaluate_basis)s
}

void %(classname)s::evaluate_basis_all(double * values,
                                       const double * x,
                                       const double * coordinate_dofs,
                                       int cell_orientation,
                                       const ufc::coordinate_mapping * cm
                                       ) const
{
%(evaluate_basis_all)s
}

void %(classname)s::evaluate_basis_derivatives(std::size_t i,
                                               std::size_t n,
                                               double * values,
                                               const double * x,
                                               const double * coordinate_dofs,
                                               int cell_orientation,
                                               const ufc::coordinate_mapping * cm
                                               ) const
{
%(evaluate_basis_derivatives)s
}

void %(classname)s::evaluate_basis_derivatives_all(std::size_t n,
                                                   double * values,
                                                   const double * x,
                                                   const double * coordinate_dofs,
                                                   int cell_orientation,
                                                   const ufc::coordinate_mapping * cm
                                                   ) const
{
%(evaluate_basis_derivatives_all)s
}

double %(classname)s::evaluate_dof(std::size_t i,
                                   const ufc::function& f,
                                   const double * coordinate_dofs,
                                   int cell_orientation,
                                   const ufc::cell& c,
                                   const ufc::coordinate_mapping * cm
                                   ) const
{
%(evaluate_dof)s
}

void %(classname)s::evaluate_dofs(double * values,
                                  const ufc::function& f,
                                  const double * coordinate_dofs,
                                  int cell_orientation,
                                  const ufc::cell& c,
                                  const ufc::coordinate_mapping * cm
                                  ) const
{
%(evaluate_dofs)s
}

void %(classname)s::interpolate_vertex_values(double * vertex_values,
                                              const double * dof_values,
                                              const double * coordinate_dofs,
                                              int cell_orientation,
                                              const ufc::coordinate_mapping * cm
                                              ) const
{
%(interpolate_vertex_values)s
}

void %(classname)s::tabulate_dof_coordinates(double * dof_coordinates,
                                             const double * coordinate_dofs,
                                             const ufc::coordinate_mapping * cm
                                             ) const
{
%(tabulate_dof_coordinates)s
}

void %(classname)s::tabulate_reference_dof_coordinates(double * reference_dof_coordinates) const
{
%(tabulate_reference_dof_coordinates)s
}

std::size_t %(classname)s::num_sub_elements() const
{
%(num_sub_elements)s
}

ufc::finite_element * %(classname)s::create_sub_element(std::size_t i) const
{
%(create_sub_element)s
}

ufc::finite_element * %(classname)s::create() const
{
%(create)s
}
"""
