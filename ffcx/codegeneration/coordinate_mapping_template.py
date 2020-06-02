# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
ufc_coordinate_mapping* create_{factory_name}(void);

// Helper used to create coordinate map using name given to the
// UFL file.
// This helper is called in user c++ code.
//
ufc_coordinate_mapping* create_coordinate_map_{prefix}(void);
"""

# declaration = """
# ufc_coordinate_mapping* create_{factory_name}(void);
# """

factory = """
// Code for coordinate mapping {factory_name}

{evaluate_reference_basis_derivatives_declaration}
void compute_jacobians_{factory_name}(double* restrict J, int num_points,
                                      const double* restrict X,
                                      const double* restrict coordinate_dofs)
{{
{compute_jacobians}
}}

void compute_jacobian_determinants_{factory_name}(double* restrict detJ, int num_points,
                                                  const double* restrict J)
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
                                     const double* restrict coordinate_dofs)
{{
{compute_geometry}
}}

void compute_midpoint_geometry_{factory_name}(double* restrict x, double* restrict J,
                                              const double* restrict coordinate_dofs)
{{
{compute_midpoint_geometry}
}}

int compute_reference_coordinates_{factory_name}(double* restrict X, int num_points,
                                                  const double* restrict x,
                                                  const double* restrict coordinate_dofs)
{{
{compute_reference_coordinates}
}}

int compute_reference_geometry_{factory_name}(double* restrict X, double* restrict J,
                                               double* restrict detJ, double* restrict K,
                                               int num_points, const double* restrict x,
                                               const double* restrict coordinate_dofs)
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
  cmap->create_scalar_dofmap = create_{scalar_dofmap_name};
  cmap->compute_physical_coordinates = compute_physical_coordinates_{factory_name};
  cmap->compute_reference_coordinates = compute_reference_coordinates_{factory_name};
  cmap->compute_reference_geometry = compute_reference_geometry_{factory_name};
  cmap->compute_jacobians = compute_jacobians_{factory_name};
  cmap->compute_jacobian_determinants = compute_jacobian_determinants_{factory_name};
  cmap->compute_jacobian_inverses = compute_jacobian_inverses_{factory_name};
  cmap->compute_geometry = compute_geometry_{factory_name};
  cmap->compute_midpoint_geometry = compute_midpoint_geometry_{factory_name};
  cmap->evaluate_reference_basis = evaluate_reference_basis_{coord_element_factory_name};
  cmap->evaluate_reference_basis_derivatives = evaluate_reference_basis_derivatives_{coord_element_factory_name};
  cmap->reference_midpoint[0] = {reference_midpoint[0]};
  cmap->reference_midpoint[1] = {reference_midpoint[1]};
  cmap->reference_midpoint[2] = {reference_midpoint[2]};
  return cmap;
}}

ufc_coordinate_mapping* create_coordinate_map_{prefix}(void)
{{
  return create_{factory_name}();
}}


// End of code for coordinate mapping {factory_name}
"""
