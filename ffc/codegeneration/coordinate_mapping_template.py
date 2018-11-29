# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
ufc_coordinate_mapping* create_{factory_name}(void);
"""

factory = """
// Code for coordinate mapping {factory_name}

{coordinate_finite_element_declaration}
ufc_finite_element* create_coordinate_finite_element_{factory_name}(void)
{{
{create_coordinate_finite_element}
}}

{coordinate_dofmap_declaration}
ufc_dofmap* create_coordinate_dofmap_{factory_name}(void)
{{
{create_coordinate_dofmap}
}}

{evaluate_reference_basis_derivatives_declaration}
void compute_jacobians_{factory_name}(double* restrict J, int num_points,
                                      const double* restrict X,
                                      const double* restrict coordinate_dofs)
{{
{compute_jacobians}
}}

void compute_jacobian_determinants_{factory_name}(double* restrict detJ, int num_points,
                                                  const double* restrict J, int cell_orientation)
{{
{compute_jacobian_determinants}
}}

void compute_jacobian_inverses_{factory_name}(double* restrict K, int num_points,
    const double* restrict J, const double* restrict detJ)
{{
{compute_jacobian_inverses}
}}

{evaluate_reference_basis_declaration}
void compute_physical_coordinates_{factory_name}(double* restrict x, int num_points,
                                                 const double* restrict X,
                                                 const double* restrict coordinate_dofs)
{{
{compute_physical_coordinates}
}}

void compute_geometry_{factory_name}(double* restrict x, double* restrict J,
                                     double* restrict detJ, double* restrict K,
                                     int num_points, const double* restrict X,
                                     const double* restrict coordinate_dofs,
                                     int cell_orientation)
{{
{compute_geometry}
}}

void compute_midpoint_geometry_{factory_name}(double* restrict x, double* restrict J,
                                              const double* restrict coordinate_dofs)
{{
{compute_midpoint_geometry}
}}

void compute_reference_coordinates_{factory_name}(double* restrict X, int num_points,
                                                  const double* restrict x,
                                                  const double* restrict coordinate_dofs,
                                                  int cell_orientation)
{{
{compute_reference_coordinates}
}}

void compute_reference_geometry_{factory_name}(double* restrict X, double* restrict J,
                                               double* restrict detJ, double* restrict K,
                                               int num_points, const double* restrict x,
                                               const double* restrict coordinate_dofs,
                                               int cell_orientation)
{{
{compute_reference_geometry}
}}




ufc_coordinate_mapping* create_{factory_name}(void)
{{
  ufc_coordinate_mapping* cmap = malloc(sizeof(*cmap));
  cmap->signature = {signature};
  cmap->create = create_{factory_name};
  cmap->geometric_dimension = {geometric_dimension};
  cmap->topological_dimension = {topological_dimension};
  cmap->cell_shape = {cell_shape};
  cmap->create_coordinate_finite_element = create_coordinate_finite_element_{factory_name};
  cmap->create_coordinate_dofmap = create_coordinate_dofmap_{factory_name};
  cmap->compute_physical_coordinates = compute_physical_coordinates_{factory_name};
  cmap->compute_reference_coordinates = compute_reference_coordinates_{factory_name};
  cmap->compute_reference_geometry = compute_reference_geometry_{factory_name};
  cmap->compute_jacobians = compute_jacobians_{factory_name};
  cmap->compute_jacobian_determinants = compute_jacobian_determinants_{factory_name};
  cmap->compute_jacobian_inverses = compute_jacobian_inverses_{factory_name};
  cmap->compute_geometry = compute_geometry_{factory_name};
  cmap->compute_midpoint_geometry = compute_midpoint_geometry_{factory_name};
  return cmap;
}}

// End of code for coordinate mapping {factory_name}
"""
