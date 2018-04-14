# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

ufc_coordinate_mapping_declaration = """
ufc_coordinate_mapping* create_{factory_name}();
"""

ufc_coordinate_mapping_factory = """
// Code for coordinate mapping {factory_name}

ufc_finite_element* create_coordinate_finite_element_{factory_name}()
{{
{create_coordinate_finite_element}
}}

ufc_dofmap* create_coordinate_dofmap_{factory_name}()
{{
{create_coordinate_dofmap}
}}

void compute_jacobians_{factory_name}(double * J, int64_t num_points,
                                      const double * X,
                                      const double * coordinate_dofs)
{{
{compute_jacobians}
}}

void compute_jacobian_determinants_{factory_name}(double * detJ, int64_t num_points,
                                                  const double * J, int cell_orientation)
{{
{compute_jacobian_determinants}
}}

void compute_jacobian_inverses_{factory_name}(
    double * K, int64_t num_points,
    const double * J, const double * detJ)
{{
{compute_jacobian_inverses}
}}

void compute_physical_coordinates_{factory_name}(double * x, int64_t num_points,
                                                 const double * X, const double * coordinate_dofs)
{{
{compute_physical_coordinates}
}}

void compute_reference_coordinates_{factory_name}(
    double * X, int64_t num_points,
    const double * x,
    const double * coordinate_dofs, int cell_orientation)
{{
{compute_reference_coordinates}
}}

void compute_reference_geometry_{factory_name}(
    double * X, double * J, double * detJ, double * K, int64_t num_points,
    const double * x,
    const double * coordinate_dofs, int cell_orientation)
{{
{compute_reference_geometry}
}}

void compute_geometry_{factory_name}(
    double * x, double * J, double * detJ, double * K, int64_t num_points,
    const double * X,
    const double * coordinate_dofs, int cell_orientation)
{{
{compute_geometry}
}}

void compute_midpoint_geometry_{factory_name}(double * x, double * J,
                                              const double * coordinate_dofs)
{{
{compute_midpoint_geometry}
}}


ufc_coordinate_mapping* create_{factory_name}()
{{
  ufc_coordinate_mapping* cmap = malloc(sizeof(*cmap));
  const char* signature = {signature};
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
