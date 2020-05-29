# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2020.

declaration = """
int init_{factory_name}(ufc_form* form);
void destroy_{factory_name}(ufc_form* form);
ufc_form* create_{factory_name}(void);

// Helper functions used to create forms based on
// the name which was given to the form in the UFL file.
// These helpers are called in user c++ code.
//
int init_{name_from_uflfile}(ufc_form* form);
void destroy_{name_from_uflfile}(ufc_form* form);
ufc_form* create_{name_from_uflfile}(void);

// Helper used to create function space using its name, as given in
// the UFL file.
// Note: There are in general more function spaces associated to the form,
//       which is the reason why the helpers takes name as argument.
// This helper is called in user c++ code.
//
ufc_function_space* create_functionspace_{name_from_uflfile}(const char* fs_name);

"""

factory = """
// Code for form {factory_name}

int original_coefficient_position_{factory_name}(int i)
{{
{original_coefficient_position}
}}

// Return a list of the coefficient names.
const char** coefficient_name_{factory_name}(void)
{{
{coefficient_name_map}
}}

// Return a list of the constant names.
const char** constant_name_{factory_name}(void)
{{
{constant_name_map}
}}

{coordinate_mapping_declaration}
int init_coordinate_mapping_{factory_name}(ufc_coordinate_mapping* cmap)
{{
  {init_coordinate_mapping}
}}

ufc_coordinate_mapping* create_coordinate_mapping_{factory_name}(void)
{{
{create_coordinate_mapping}
}}

{finite_element_declaration}
int init_finite_element_{factory_name}(ufc_finite_element* element, int i)
{{
  {init_finite_element}
}}

ufc_finite_element* create_finite_element_{factory_name}(int i)
{{
{create_finite_element}
}}

{dofmap_declaration}
int init_dofmap_{factory_name}(ufc_dofmap* dofmap, int i)
{{
  {init_dofmap}
}}

ufc_dofmap* create_dofmap_{factory_name}(int i)
{{
{create_dofmap}
}}

int init_cell_integral_{factory_name}(ufc_integral* integral, int subdomain_id)
{{
  {init_cell_integral}
}}

ufc_integral* create_cell_integral_{factory_name}(int subdomain_id)
{{
  {create_cell_integral}
}}

void get_cell_integral_ids_{factory_name}(int *ids)
{{
  {get_cell_integral_ids}
}}

int init_exterior_facet_integral_{factory_name}(ufc_integral* integral, int subdomain_id)
{{
  {init_exterior_facet_integral}
}}

ufc_integral* create_exterior_facet_integral_{factory_name}(int subdomain_id)
{{
  {create_exterior_facet_integral}
}}

void get_exterior_facet_integral_ids_{factory_name}(int *ids)
{{
  {get_exterior_facet_integral_ids}
}}

int init_interior_facet_integral_{factory_name}(ufc_integral* integral, int subdomain_id)
{{
  {init_interior_facet_integral}
}}

ufc_integral* create_interior_facet_integral_{factory_name}(int subdomain_id)
{{
{create_interior_facet_integral}
}}

void get_interior_facet_integral_ids_{factory_name}(int *ids)
{{
  {get_interior_facet_integral_ids}
}}

int init_vertex_integral_{factory_name}(ufc_integral* integral, int subdomain_id)
{{
  {init_vertex_integral}
}}

ufc_integral* create_vertex_integral_{factory_name}(int subdomain_id)
{{
{create_vertex_integral}
}}

void get_vertex_integral_ids_{factory_name}(int *ids)
{{
  {get_vertex_integral_ids}
}}

int init_custom_integral_{factory_name}(ufc_custom_integral* integral, int subdomain_id)
{{
  {init_custom_integral}
}}

ufc_custom_integral* create_custom_integral_{factory_name}(int subdomain_id)
{{
{create_custom_integral}
}}

void get_custom_integral_ids_{factory_name}(int *ids)
{{
  {get_custom_integral_ids}
}}

void destroy_{factory_name}(ufc_form* form);
ufc_form* create_{factory_name}(void);

int init_{factory_name}(ufc_form* form)
{{
  form->signature = {signature};
  form->rank = {rank};
  form->num_coefficients = {num_coefficients};
  form->num_constants = {num_constants};

  form->init = init_{factory_name};
  form->destroy = destroy_{factory_name};
  form->create = create_{factory_name};

  form->original_coefficient_position = original_coefficient_position_{factory_name};

  form->coefficient_name_map = coefficient_name_{factory_name};
  form->constant_name_map = constant_name_{factory_name};

  form->init_coordinate_mapping = init_coordinate_mapping_{factory_name};
  form->create_coordinate_mapping = create_coordinate_mapping_{factory_name};
  form->init_finite_element = init_finite_element_{factory_name};
  form->create_finite_element = create_finite_element_{factory_name};
  form->init_dofmap = init_dofmap_{factory_name};
  form->create_dofmap = create_dofmap_{factory_name};

  form->get_cell_integral_ids = get_cell_integral_ids_{factory_name};
  form->get_exterior_facet_integral_ids = get_exterior_facet_integral_ids_{factory_name};
  form->get_interior_facet_integral_ids = get_interior_facet_integral_ids_{factory_name};
  form->get_vertex_integral_ids = get_vertex_integral_ids_{factory_name};
  form->get_custom_integral_ids = get_custom_integral_ids_{factory_name};

  form->num_cell_integrals = {num_cell_integrals};
  form->num_exterior_facet_integrals = {num_exterior_facet_integrals};
  form->num_interior_facet_integrals = {num_interior_facet_integrals};
  form->num_vertex_integrals = {num_vertex_integrals};
  form->num_custom_integrals = {num_custom_integrals};

  form->init_cell_integral = init_cell_integral_{factory_name};
  form->create_cell_integral = create_cell_integral_{factory_name};
  form->init_exterior_facet_integral = init_exterior_facet_integral_{factory_name};
  form->create_exterior_facet_integral = create_exterior_facet_integral_{factory_name};
  form->init_interior_facet_integral = init_interior_facet_integral_{factory_name};
  form->create_interior_facet_integral = create_interior_facet_integral_{factory_name};
  form->init_vertex_integral = init_vertex_integral_{factory_name};
  form->create_vertex_integral = create_vertex_integral_{factory_name};
  form->init_custom_integral = init_custom_integral_{factory_name};
  form->create_custom_integral = create_custom_integral_{factory_name};
  return 0;
}}

void destroy_{factory_name}(ufc_form* form)
{{
}}

ufc_form* create_{factory_name}(void)
{{
  ufc_form* form = malloc(sizeof(*form));
  init_{factory_name}(form);
  return form;
}}

int init_{name_from_uflfile}(ufc_form* form)
{{
  return init_{factory_name}(form);
}}

void destroy_{name_from_uflfile}(ufc_form* form)
{{
  destroy_{factory_name}(form);
}}

ufc_form* create_{name_from_uflfile}(void)
{{
  return create_{factory_name}();
}}

ufc_function_space* create_functionspace_{name_from_uflfile}(const char* function_name)
{{
  {create_functionspace}
}}

// End of code for form {factory_name}
"""
