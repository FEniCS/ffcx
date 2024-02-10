# Copyright (C) 2020 Matthew Scroggs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import basix.ufl
import numpy as np

import ffcx
import ffcx.codegeneration.jit


def test_finite_element(compile_args):
    ufl_element = basix.ufl.element("Lagrange", "triangle", 1)
    jit_compiled_elements, module, code = ffcx.codegeneration.jit.compile_elements(
        [ufl_element], cffi_extra_compile_args=compile_args
    )
    ufcx_element, ufcx_dofmap = jit_compiled_elements[0]

    assert ufcx_element.topological_dimension == 2
    assert ufcx_element.space_dimension == 3
    assert ufcx_element.reference_value_rank == 0
    assert ufcx_element.reference_value_size == 1
    assert ufcx_element.block_size == 1
    assert ufcx_element.num_sub_elements == 0

    assert ufcx_dofmap.block_size == 1
    assert ufcx_dofmap.num_global_support_dofs == 0
    assert ufcx_dofmap.num_global_support_dofs == 0
    assert ufcx_dofmap.num_element_support_dofs == 3
    off = np.array([ufcx_dofmap.entity_dof_offsets[i] for i in range(8)])
    assert np.all(np.diff(off) == [1, 1, 1, 0, 0, 0, 0])

    for v in range(3):
        assert ufcx_dofmap.entity_dofs[v] == v
    assert ufcx_dofmap.num_sub_dofmaps == 0


def test_vector_element(compile_args):
    ufl_element = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
    jit_compiled_elements, module, code = ffcx.codegeneration.jit.compile_elements(
        [ufl_element], cffi_extra_compile_args=compile_args
    )
    ufcx_element, ufcx_dofmap = jit_compiled_elements[0]

    assert ufcx_element.topological_dimension == 2
    assert ufcx_element.space_dimension == 6
    assert ufcx_element.reference_value_rank == 1
    assert ufcx_element.reference_value_shape[0] == 2
    assert ufcx_element.reference_value_size == 2
    assert ufcx_element.block_size == 2
    assert ufcx_element.num_sub_elements == 2

    assert ufcx_dofmap.block_size == 2
    assert ufcx_dofmap.num_global_support_dofs == 0
    assert ufcx_dofmap.num_global_support_dofs == 0
    assert ufcx_dofmap.num_element_support_dofs == 3
    off = np.array([ufcx_dofmap.entity_dof_offsets[i] for i in range(8)])
    assert np.all(np.diff(off) == [1, 1, 1, 0, 0, 0, 0])

    for v in range(3):
        assert ufcx_dofmap.entity_dofs[v] == v
    assert ufcx_dofmap.num_sub_dofmaps == 2


def test_tensor_element(compile_args):
    ufl_element = basix.ufl.element("Lagrange", "triangle", 1, shape=(2, 2))
    jit_compiled_elements, module, code = ffcx.codegeneration.jit.compile_elements(
        [ufl_element], cffi_extra_compile_args=compile_args
    )
    ufcx_element, ufcx_dofmap = jit_compiled_elements[0]

    assert ufcx_element.topological_dimension == 2
    assert ufcx_element.space_dimension == 12
    assert ufcx_element.reference_value_rank == 2
    assert ufcx_element.reference_value_shape[0] == 2
    assert ufcx_element.reference_value_shape[1] == 2
    assert ufcx_element.reference_value_size == 4
    assert ufcx_element.block_size == 4
    assert ufcx_element.num_sub_elements == 4

    assert ufcx_dofmap.block_size == 4
    assert ufcx_dofmap.num_global_support_dofs == 0
    assert ufcx_dofmap.num_global_support_dofs == 0
    assert ufcx_dofmap.num_element_support_dofs == 3
    off = np.array([ufcx_dofmap.entity_dof_offsets[i] for i in range(8)])
    assert np.all(np.diff(off) == [1, 1, 1, 0, 0, 0, 0])

    for v in range(3):
        assert ufcx_dofmap.entity_dofs[v] == v
    assert ufcx_dofmap.num_sub_dofmaps == 4


def test_vector_quadrature_element(compile_args):
    ufl_element = basix.ufl.blocked_element(
        basix.ufl.quadrature_element("tetrahedron", degree=2, scheme="default"), shape=(3,)
    )
    jit_compiled_elements, module, code = ffcx.codegeneration.jit.compile_elements(
        [ufl_element], cffi_extra_compile_args=compile_args
    )
    ufcx_element, ufcx_dofmap = jit_compiled_elements[0]

    assert ufcx_element.topological_dimension == 3
    assert ufcx_element.space_dimension == 12
    assert ufcx_element.reference_value_rank == 1
    assert ufcx_element.reference_value_shape[0] == 3
    assert ufcx_element.reference_value_size == 3
    assert ufcx_element.block_size == 3
    assert ufcx_element.num_sub_elements == 3

    assert ufcx_dofmap.block_size == 3
    assert ufcx_dofmap.num_global_support_dofs == 0
    assert ufcx_dofmap.num_global_support_dofs == 0
    assert ufcx_dofmap.num_element_support_dofs == 4
    off = np.array([ufcx_dofmap.entity_dof_offsets[i] for i in range(16)])
    assert np.all(np.diff(off) == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4])

    for i in range(4):
        assert ufcx_dofmap.entity_dofs[i] == i

    assert ufcx_dofmap.num_sub_dofmaps == 3
