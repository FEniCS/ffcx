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

UFC_DECL = """
typedef enum
{
interval = 10,
triangle = 20,
quadrilateral = 30,
tetrahedron = 40,
hexahedron = 50,
vertex = 60,
} ufc_shape;

typedef struct ufc_coordinate_mapping ufc_coordinate_mapping;
typedef struct ufc_finite_element ufc_finite_element;

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
void (*transform_values)(
    double* restrict reference_values,
    const double* restrict physical_values,
    const double* restrict coordinate_dofs,
    int cell_orientation, const ufc_coordinate_mapping* cm);
void (*tabulate_reference_dof_coordinates)(
    double* restrict reference_dof_coordinates);
int num_sub_elements;
ufc_finite_element* (*create_sub_element)(int i);
ufc_finite_element* (*create)(void);
} ufc_finite_element;

"""


def compile_elements(elements, module_name=None):
    """Compile a list of UFL elements into UFC Python objects"""
    code_body = ""
    decl = UFC_DECL
    element_template = "ufc_finite_element * create_{name}(void);"
    for e in elements:
        _, impl = ffc.compiler.compile_element(e)
        code_body += impl.replace("#include \"Element.h\"", "#include <ufc.h>")

        p = ffc.parameters.validate_parameters(None)
        name = ffc.representation.make_finite_element_jit_classname(e, p)
        create_element = element_template.format(name=name)
        decl += create_element + "\n"

    if not module_name:
        h = hashlib.sha1()
        h.update((code_body + decl).encode('utf-8'))
        module_name = "_" + h.hexdigest()

    ffibuilder = cffi.FFI()
    ffibuilder.set_source(
        module_name, code_body, include_dirs=[ffc.backends.ufc.get_include_path()])
    ffibuilder.cdef(decl)

    compile_dir = "compile_cache"
    ffibuilder.compile(tmpdir=compile_dir, verbose=False)

    # Build list of compiled elements
    compiled_elements = []
    compiled_module = importlib.import_module(compile_dir + "." + module_name)
    for e in elements:
        p = ffc.parameters.validate_parameters(None)
        name = ffc.representation.make_finite_element_jit_classname(e, p)
        create_element = "create_" + name
        compiled_elements.append(getattr(compiled_module.lib, create_element)())

    return compiled_elements, compiled_module
