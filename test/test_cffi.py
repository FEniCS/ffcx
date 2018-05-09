# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
from cffi import FFI

import ffc
import ufl

ufc_decl = """
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
void (*map_dofs)(double* restrict values, const double* restrict vals,
                    const double* restrict coordinate_dofs,
                    int cell_orientation, const ufc_coordinate_mapping* cm);
void (*tabulate_reference_dof_coordinates)(
    double* restrict reference_dof_coordinates);
int num_sub_elements;
ufc_finite_element* (*create_sub_element)(int i);
ufc_finite_element* (*create)(void);
} ufc_finite_element;

"""

element_template = "ufc_finite_element * create_{name}(void);"


def jit_elements(elements):
    code = ""
    decl = ufc_decl
    for e in elements:
        header, impl = ffc.compiler.compile_element(e)
        code += impl.replace("#include \"Element.h\"", "#include <ufc.h>")

        p = ffc.parameters.validate_parameters(None)
        name = ffc.representation.make_finite_element_jit_classname(e, p)
        create_element = element_template.format(name=name)
        decl += create_element + "\n"

    ufc_path = ffc.backends.ufc.get_include_path()
    ffibuilder = FFI()
    ffibuilder.set_source("_test", code, include_dirs=[ufc_path])
    ffibuilder.cdef(decl)
    ffibuilder.compile(verbose=True)

    compiled_elements = []
    import _test
    from _test import ffi, lib
    print(ffi)
    print(lib)
    for e in elements:
        p = ffc.parameters.validate_parameters(None)
        name = ffc.representation.make_finite_element_jit_classname(e, p)
        create_element = "create_" + name
        print(create_element)
        compiled_elements.append(getattr(_test.lib, create_element)())

    return compiled_elements, _test


# Build list of UFL elements
cell = ufl.triangle
elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 3)]

# Compile elements
compiled_elements, module = jit_elements(elements)


# Test
for e, compiled_e in zip(elements, compiled_elements):
    print(compiled_e.geometric_dimension)
    print(compiled_e.degree)
    print(compiled_e)
    print(compiled_e.family)
    test = module.ffi.string(compiled_e.family)

    tdim = compiled_e.topological_dimension
    space_dim = compiled_e.space_dimension
    X = np.zeros([space_dim, tdim])

    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_e.tabulate_reference_dof_coordinates(X_ptr)
    print(X)
