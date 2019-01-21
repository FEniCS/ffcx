# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import hashlib
import importlib

import cffi

import ffc

UFC_HEADER_DECL = """
typedef double {} ufc_scalar_t;  /* Hack to deal with scalar type */

typedef struct ufc_coordinate_mapping ufc_coordinate_mapping;
typedef struct ufc_finite_element ufc_finite_element;
typedef struct ufc_dofmap ufc_dofmap;

typedef enum
{{
interval = 10,
triangle = 20,
quadrilateral = 30,
tetrahedron = 40,
hexahedron = 50,
vertex = 60,
}} ufc_shape;
"""

UFC_ELEMENT_DECL = """
typedef struct ufc_finite_element
{
const char* signature;
ufc_shape cell_shape;
int topological_dimension;
int geometric_dimension;
int space_dimension;
int value_rank;
int (*value_dimension)(int i);
int value_size;
int reference_value_rank;
int (*reference_value_dimension)(int i);
int reference_value_size;
int degree;
const char* family;
int (*evaluate_reference_basis)(double* restrict reference_values,
                                int num_points, const double* restrict X);
int (*evaluate_reference_basis_derivatives)(
    double* restrict reference_values, int order, int num_points,
    const double* restrict X);
int (*transform_reference_basis_derivatives)(
    double* restrict values, int order, int num_points,
    const double* restrict reference_values, const double* restrict X,
    const double* restrict J, const double* restrict detJ,
    const double* restrict K, int cell_orientation);
int (*transform_values)(
    ufc_scalar_t* restrict reference_values,
    const ufc_scalar_t* restrict physical_values,
    const double* restrict coordinate_dofs,
    int cell_orientation, const ufc_coordinate_mapping* cm);
int (*tabulate_reference_dof_coordinates)(
    double* restrict reference_dof_coordinates);
int num_sub_elements;
ufc_finite_element* (*create_sub_element)(int i);
ufc_finite_element* (*create)(void);
} ufc_finite_element;
"""

UFC_DOFMAP_DECL = """
typedef struct ufc_dofmap
{
const char* signature;
int num_global_support_dofs;
int num_element_support_dofs;
int num_entity_dofs[4];
int num_entity_closure_dofs[4];
void (*tabulate_dof_permutations)(int* restrict perm, const int64_t* restrict global_indices);
void (*tabulate_entity_dofs)(int* restrict dofs, int d, int i);
void (*tabulate_entity_closure_dofs)(int* restrict dofs, int d, int i);
int num_sub_dofmaps;
ufc_dofmap* (*create_sub_dofmap)(int i);
ufc_dofmap* (*create)(void);
} ufc_dofmap;
"""

UFC_COORDINATEMAPPING_DECL = """
typedef struct ufc_coordinate_mapping
{
const char* signature;
ufc_coordinate_mapping* (*create)(void);
int geometric_dimension;
int topological_dimension;
ufc_shape cell_shape;
ufc_finite_element* (*create_coordinate_finite_element)(void);
ufc_dofmap* (*create_coordinate_dofmap)(void);
void (*compute_physical_coordinates)(
    double* restrict x, int num_points, const double* restrict X,
    const double* restrict coordinate_dofs);
void (*compute_reference_coordinates)(
    double* restrict X, int num_points, const double* restrict x,
    const double* restrict coordinate_dofs, int cell_orientation);
void (*compute_reference_geometry)(double* restrict X, double* restrict J,
                                    double* restrict detJ,
                                    double* restrict K, int num_points,
                                    const double* restrict x,
                                    const double* restrict coordinate_dofs,
                                    int cell_orientation);
void (*compute_jacobians)(double* restrict J, int num_points,
                            const double* restrict X,
                            const double* restrict coordinate_dofs);
void (*compute_jacobian_determinants)(double* restrict detJ, int num_points,
                                        const double* restrict J,
                                        int cell_orientation);
void (*compute_jacobian_inverses)(double* restrict K, int num_points,
                                    const double* restrict J,
                                    const double* restrict detJ);
void (*compute_geometry)(double* restrict x, double* restrict J,
                            double* restrict detJ, double* restrict K,
                            int num_points, const double* restrict X,
                            const double* restrict coordinate_dofs,
                            int cell_orientation);
void (*compute_midpoint_geometry)(double* restrict x, double* restrict J,
                                    const double* restrict coordinate_dofs);

} ufc_coordinate_mapping;
"""

UFC_INTEGRAL_DECL = """
typedef struct ufc_cell_integral
{
const bool* enabled_coefficients;
void (*tabulate_tensor)(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                        const double* restrict coordinate_dofs,
                        int cell_orientation);
} ufc_cell_integral;

typedef struct ufc_exterior_facet_integral
{
const bool* enabled_coefficients;
void (*tabulate_tensor)(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                        const double* restrict coordinate_dofs, int facet,
                        int cell_orientation);
} ufc_exterior_facet_integral;

typedef struct ufc_interior_facet_integral
{
const bool* enabled_coefficients;
void (*tabulate_tensor)(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                        const double* restrict coordinate_dofs_0,
                        const double* restrict coordinate_dofs_1,
                        int facet_0, int facet_1, int cell_orientation_0,
                        int cell_orientation_1);
} ufc_interior_facet_integral;

typedef struct ufc_vertex_integral
{
const bool* enabled_coefficients;
void (*tabulate_tensor)(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                        const double* restrict coordinate_dofs, int vertex,
                        int cell_orientation);
} ufc_vertex_integral;

typedef struct ufc_custom_integral
{
const bool* enabled_coefficients;
void (*tabulate_tensor)(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                        const double* restrict coordinate_dofs,
                        int num_quadrature_points,
                        const double* restrict quadrature_points,
                        const double* restrict quadrature_weights,
                        const double* restrict facet_normals,
                        int cell_orientation);
} ufc_custom_integral;
"""

UFC_FORM_DECL = """
typedef struct ufc_form
{
const char* signature;
int rank;
int num_coefficients;
int (*original_coefficient_position)(int i);
ufc_finite_element* (*create_coordinate_finite_element)(void);
ufc_dofmap* (*create_coordinate_dofmap)(void);
ufc_coordinate_mapping* (*create_coordinate_mapping)(void);
ufc_finite_element* (*create_finite_element)(int i);
ufc_dofmap* (*create_dofmap)(int i);
int max_cell_subdomain_id;
int max_exterior_facet_subdomain_id;
int max_interior_facet_subdomain_id;
int max_vertex_subdomain_id;
int max_custom_subdomain_id;
bool has_cell_integrals;
bool has_exterior_facet_integrals;
bool has_interior_facet_integrals;
bool has_vertex_integrals;
bool has_custom_integrals;
ufc_cell_integral* (*create_cell_integral)(int subdomain_id);
ufc_exterior_facet_integral* (*create_exterior_facet_integral)(
    int subdomain_id);
ufc_interior_facet_integral* (*create_interior_facet_integral)(
    int subdomain_id);
ufc_vertex_integral* (*create_vertex_integral)(int subdomain_id);
ufc_custom_integral* (*create_custom_integral)(int subdomain_id);
ufc_cell_integral* (*create_default_cell_integral)(void);
ufc_exterior_facet_integral* (*create_default_exterior_facet_integral)(void);
ufc_interior_facet_integral* (*create_default_interior_facet_integral)(void);
ufc_vertex_integral* (*create_default_vertex_integral)(void);
ufc_custom_integral* (*create_default_custom_integral)(void);
} ufc_form;
"""


def compile_elements(elements, module_name=None, parameters=None):
    """Compile a list of UFL elements and dofmaps into UFC Python objects"""
    p = ffc.parameters.validate_parameters(parameters)

    decl = UFC_HEADER_DECL.format("") + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL
    element_template = "ufc_finite_element * create_{name}(void);\n"
    dofmap_template = "ufc_dofmap * create_{name}(void);\n"
    names = []
    for e in elements:
        name = ffc.representation.make_finite_element_jit_classname(e, p)
        names.append(name)
        decl += element_template.format(name=name)
#        name = ffc.representation.make_dofmap_jit_classname(e, p)
#        names.append(name)
#        decl += dofmap_template.format(name=name)

    _, code_body = ffc.compiler.compile_element(elements, parameters=p)

    return _compile_objects(decl, code_body, names, module_name, p)


def compile_forms(forms, module_name=None, parameters=None):
    """Compile a list of UFL forms into UFC Python objects"""
    p = ffc.parameters.validate_parameters(parameters)

    if p and "complex" in p["scalar_type"]:
        complex_mode = "_Complex"
    else:
        complex_mode = ""
    decl = UFC_HEADER_DECL.format(complex_mode) + UFC_ELEMENT_DECL \
        + UFC_DOFMAP_DECL + UFC_COORDINATEMAPPING_DECL \
        + UFC_INTEGRAL_DECL + UFC_FORM_DECL

    form_names = [ffc.classname.make_name("Form", "form", i)
                  for i in range(len(forms))]
    form_template = "ufc_form * create_{name}(void);\n"
    for name in form_names:
        decl += form_template.format(name=name)

    _, code_body = ffc.compiler.compile_form(forms, parameters=p)

    return _compile_objects(decl, code_body, form_names, module_name, p)


def compile_coordinate_maps(cmaps, module_name=None, parameters=None):
    """Compile a list of UFL coordinate mappings into UFC Python objects"""
    p = ffc.parameters.validate_parameters(parameters)

    decl = UFC_HEADER_DECL.format("") + UFC_COORDINATEMAPPING_DECL
    cmap_template = "ufc_coordinate_mapping * create_{name}(void);\n"
    cmap_names = [ffc.representation.make_coordinate_mapping_jit_classname(cmap, p) for cmap in cmaps]
    for name in cmap_names:
        decl += cmap_template.format(name=name)

    _, code_body = ffc.compiler.compile_coordinate_mapping(cmaps, parameters=p)

    return _compile_objects(decl, code_body, cmap_names, module_name, p)


def _compile_objects(decl, code_body, object_names, module_name, parameters):

    if not module_name:
        h = hashlib.sha1()
        h.update((code_body + decl).encode('utf-8'))
        module_name = "_" + h.hexdigest()

    ffibuilder = cffi.FFI()
    ffibuilder.set_source(
        module_name, code_body, include_dirs=[ffc.codegeneration.get_include_path()])
    ffibuilder.cdef(decl)

    cache_dir = None
    if parameters:
        cache_dir = parameters.get("cache_dir")
    if not cache_dir:
        cache_dir = "compile_cache"

    ffibuilder.compile(tmpdir=cache_dir, verbose=False)

    # Build list of compiled objects
    compiled_module = importlib.import_module(cache_dir + "." + module_name)
    print(object_names)
    compiled_objects = [getattr(compiled_module.lib, "create_" + name)() for name in object_names]

    return compiled_objects, compiled_module
