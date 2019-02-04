# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import importlib
import sys
import os
import logging
import time
import cffi
import pathlib

import ufl
import ffc
from ffc.analysis import analyze_ufl_objects

logger = logging.getLogger(__name__)

UFC_HEADER_DECL = """
typedef {} ufc_scalar_t;  /* Hack to deal with scalar type */

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


def get_ufl_dependencies(ufl_objects, parameters):

    _, unique_elements, _, unique_coordinate_elements = analyze_ufl_objects(ufl_objects, parameters)

    mesh_id = None
    if isinstance(ufl_objects[0], ufl.Form):
        mesh_id = ufl_objects[0].ufl_domain().ufl_id()
    elif isinstance(ufl_objects[0], ufl.Mesh):
        mesh_id = ufl_objects[0].ufl_id()
        unique_meshes = []
    if mesh_id is not None:
        unique_meshes = [ufl.Mesh(element, ufl_id=mesh_id) for element in unique_coordinate_elements]

    # Avoid returning self as dependency for infinite recursion
    unique_elements = tuple(
        element for element in unique_elements if element not in ufl_objects)
    unique_meshes = tuple(mesh for mesh in unique_meshes if mesh not in ufl_objects)

    # Setup dependencies (these will be jitted before continuing to
    # compile ufl_objects)
    dependent_ufl_objects = {
        "element": unique_elements,
        "coordinate_mapping": unique_meshes,
    }
    return dependent_ufl_objects


def get_cached_module(module_name, object_names, parameters):
    cache_dir = pathlib.Path(parameters.get("cache_dir",
                                            "compile_cache"))
    cache_dir = cache_dir.expanduser()

    timeout = int(parameters.get("timeout", 10))

    c_filename = cache_dir.joinpath(module_name + ".c")
    ready_name = c_filename.with_suffix(".c.cached")

    # Ensure cache dir exists
    os.makedirs(cache_dir, exist_ok=True)

    # Ensure it is first on the path for loading modules
    sys.path.insert(0, str(cache_dir))

    try:
        # Create C file with exclusive access
        open(c_filename, "x")
        return None, None

    except FileExistsError:
        print("Cached C file already exists:" + str(c_filename))
        # Now wait for ready
        for i in range(timeout):
            if os.path.exists(ready_name):
                # Build list of compiled objects
                compiled_module = \
                    importlib.import_module(module_name)
                sys.path.remove(str(cache_dir))
                compiled_objects = \
                    [getattr(compiled_module.lib,
                             "create_" + name)()
                     for name in object_names]

                return compiled_objects, compiled_module

            print("Waiting for " + str(ready_name) + " to appear.")
            time.sleep(1)
        raise TimeoutError("""JIT compilation did not complete on another process.
        Try cleaning cache or increase timeout parameter.""")


def compile_elements(elements, module_name=None, parameters=None):
    """Compile a list of UFL elements and dofmaps into UFC Python objects"""
    p = ffc.parameters.validate_parameters(parameters)

    print('Compiliing element with ' + str(elements) + 'params = ' + str(p))

    # Get a signature for these elements
    module_name = 'elements_' + ffc.classname.compute_signature(elements, '', p)
    print('element module name = ', module_name)

    names = []
    for e in elements:
        name = ffc.ir.representation.make_finite_element_jit_classname(e, "JIT", p)
        names.append(name)
        name = ffc.ir.representation.make_dofmap_jit_classname(e, "JIT", p)
        names.append(name)

    obj, mod = get_cached_module(module_name, names, p)
    if obj is not None:
        print(' **** Loaded from cache: ' + module_name)
        # Pair up elements with dofmaps
        obj = list(zip(obj[::2], obj[1::2]))
        return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL
    element_template = "ufc_finite_element * create_{name}(void);\n"
    dofmap_template = "ufc_dofmap * create_{name}(void);\n"

    for i in range(len(elements)):
        decl += element_template.format(name=names[i * 2])
        decl += dofmap_template.format(name=names[i * 2 + 1])

    _, code_body = ffc.compiler.compile_ufl_objects(elements, prefix="JIT", parameters=p)

    objects, module = _compile_objects(decl, code_body, names, module_name, p)
    # Pair up elements with dofmaps
    objects = list(zip(objects[::2], objects[1::2]))
    return objects, module


def compile_forms(forms, module_name=None, parameters=None):
    """Compile a list of UFL forms into UFC Python objects"""
    p = ffc.parameters.validate_parameters(parameters)

    # Get a signature for these forms
    module_name = 'forms_' + ffc.classname.compute_signature(forms, '', p)
    print('form module name = ', module_name)

    form_names = [ffc.classname.make_name("JIT", "form", i)
                  for i in range(len(forms))]

    deps = get_ufl_dependencies(forms, p)

    print('*** DEPS = ', deps)
    for k, v in deps.items():
        if (k == 'element'):
            compile_elements(v, parameters=p)
        if (k == 'coordinate_mapping'):
            compile_coordinate_maps(v, parameters=p)

    obj, mod = get_cached_module(module_name, form_names, p)
    if obj is not None:
        print('Loaded from cache: ' + module_name)
        return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL \
        + UFC_DOFMAP_DECL + UFC_COORDINATEMAPPING_DECL \
        + UFC_INTEGRAL_DECL + UFC_FORM_DECL

    form_template = "ufc_form * create_{name}(void);\n"
    for name in form_names:
        decl += form_template.format(name=name)

    _, code_body = ffc.compiler.compile_ufl_objects(forms, prefix="JIT", parameters=p, jit=True)

    return _compile_objects(decl, code_body, form_names, module_name, p)


def compile_coordinate_maps(meshes, module_name=None, parameters=None):
    """Compile a list of UFL coordinate mappings into UFC Python objects"""
    p = ffc.parameters.validate_parameters(parameters)

    # Get a signature for these cmaps
    module_name = 'cmaps_' + ffc.classname.compute_signature(meshes, '', p, True)
    print('cmap module name = ', module_name)

    cmap_names = [ffc.ir.representation.make_coordinate_mapping_jit_classname(
        mesh.ufl_coordinate_element(), "JIT", p) for mesh in meshes]

    obj, mod = get_cached_module(module_name, cmap_names, p)
    if obj is not None:
        print('Loaded from cache: ' + module_name)
        return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_COORDINATEMAPPING_DECL
    cmap_template = "ufc_coordinate_mapping * create_{name}(void);\n"

    for name in cmap_names:
        decl += cmap_template.format(name=name)

    _, code_body = ffc.compiler.compile_ufl_objects(meshes, prefix="JIT", parameters=p)

    return _compile_objects(decl, code_body, cmap_names, module_name, p)


def _compile_objects(decl, code_body, object_names, module_name, parameters):

    ffibuilder = cffi.FFI()
    ffibuilder.set_source(
        module_name, code_body, include_dirs=[ffc.codegeneration.get_include_path()])

    ffibuilder.cdef(decl)

    cache_dir = pathlib.Path(parameters.get("cache_dir",
                                            "compile_cache"))
    cache_dir = cache_dir.expanduser()
    c_filename = cache_dir.joinpath(module_name + ".c")
    ready_name = c_filename.with_suffix(".c.cached")

    # Ensure path is set for module
    sys.path.insert(0, str(cache_dir))

    # Ensure cache dir exists
    os.makedirs(cache_dir, exist_ok=True)

    # Compile
    ffibuilder.compile(tmpdir=cache_dir, verbose=False)

    # Create a "status ready" file
    # If this fails, it is an error, because it should not exist yet.
    fd = open(ready_name, "x")
    fd.close()

    # Build list of compiled objects
    compiled_module = importlib.import_module(module_name)
    sys.path.remove(str(cache_dir))
    compiled_objects = [getattr(compiled_module.lib, "create_" + name)() for name in object_names]

    return compiled_objects, compiled_module
