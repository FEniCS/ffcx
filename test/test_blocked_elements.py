# Copyright (C) 2020 Matthew Scroggs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np

import ffcx
import ffcx.codegeneration.jit
import ufl


def test_finite_element(compile_args):
    ufl_element = ufl.FiniteElement("Lagrange", ufl.triangle, 1)
    jit_compiled_elements, module, code = ffcx.codegeneration.jit.compile_elements(
        [ufl_element], cffi_extra_compile_args=compile_args)
    ufc_element, ufc_dofmap = jit_compiled_elements[0]

    assert ufc_element.topological_dimension == 2
    assert ufc_element.geometric_dimension == 2
    assert ufc_element.space_dimension == 3
    assert ufc_element.value_rank == 0
    assert ufc_element.value_size == 1
    assert ufc_element.reference_value_rank == 0
    assert ufc_element.reference_value_size == 1
    assert ufc_element.block_size == 1
    assert ufc_element.num_sub_elements == 0

    assert ufc_dofmap.block_size == 1
    assert ufc_dofmap.num_global_support_dofs == 0
    assert ufc_dofmap.num_global_support_dofs == 0
    assert ufc_dofmap.num_element_support_dofs == 3
    assert ufc_dofmap.num_sub_dofmaps == 0


def test_vector_element(compile_args):
    ufl_element = ufl.VectorElement("Lagrange", ufl.triangle, 1)
    jit_compiled_elements, module, code = ffcx.codegeneration.jit.compile_elements(
        [ufl_element], cffi_extra_compile_args=compile_args)
    ufc_element, ufc_dofmap = jit_compiled_elements[0]

    assert ufc_element.topological_dimension == 2
    assert ufc_element.geometric_dimension == 2
    assert ufc_element.space_dimension == 6
    assert ufc_element.value_rank == 1
    assert ufc_element.value_shape[0] == 2
    assert ufc_element.value_size == 2
    assert ufc_element.reference_value_rank == 1
    assert ufc_element.reference_value_shape[0] == 2
    assert ufc_element.reference_value_size == 2
    assert ufc_element.block_size == 2
    assert ufc_element.num_sub_elements == 2

    assert ufc_dofmap.block_size == 2
    assert ufc_dofmap.num_global_support_dofs == 0
    assert ufc_dofmap.num_global_support_dofs == 0
    assert ufc_dofmap.num_element_support_dofs == 3
    assert ufc_dofmap.num_sub_dofmaps == 2


def test_tensor_element(compile_args):
    ufl_element = ufl.TensorElement("Lagrange", ufl.triangle, 1)
    jit_compiled_elements, module, code = ffcx.codegeneration.jit.compile_elements(
        [ufl_element], cffi_extra_compile_args=compile_args)
    ufc_element, ufc_dofmap = jit_compiled_elements[0]

    assert ufc_element.topological_dimension == 2
    assert ufc_element.geometric_dimension == 2
    assert ufc_element.space_dimension == 12
    assert ufc_element.value_rank == 2
    assert ufc_element.value_shape[0] == 2
    assert ufc_element.value_shape[1] == 2
    assert ufc_element.value_size == 4
    assert ufc_element.reference_value_rank == 2
    assert ufc_element.reference_value_shape[0] == 2
    assert ufc_element.reference_value_shape[1] == 2
    assert ufc_element.reference_value_size == 4
    assert ufc_element.block_size == 4
    assert ufc_element.num_sub_elements == 4

    assert ufc_dofmap.block_size == 4
    assert ufc_dofmap.num_global_support_dofs == 0
    assert ufc_dofmap.num_global_support_dofs == 0
    assert ufc_dofmap.num_element_support_dofs == 3
    assert ufc_dofmap.num_sub_dofmaps == 4
