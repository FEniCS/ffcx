# -*- coding: utf-8 -*-

import pytest

import numpy

from ffc.uflacs.backends.ufc.generators import *

import ffc.uflacs.language.cnodes as L

# from ffc.cpp import set_float_formatting
# set_float_formatting(precision=8)


# TODO: Make this a feature of dijitso: dijitso show-function modulehash functionname
def extract_function(name, code):
    lines = code.split("\n")
    n = len(lines)
    begin = None
    body = None
    for i in range(n):
        if (name + "(") in lines[i]:
            for j in range(i, n):
                if lines[j] == "{":
                    begin = i
                    body = j
                    break
            break
    if begin is None:
        return "didnt find %s" % (name,)
    end = n
    for i in range(body, n):
        if lines[i] == "}":
            end = i + 1
            break
    sublines = lines[begin:end]
    return '\n'.join(sublines)


def basic_class_properties(classname):
    ir = {
        "classname": classname,
        "constructor": "",
        "constructor_arguments": "",
        "initializer_list": "",
        "destructor": "",
        "members": "",
        "preamble": "",
        "jit": False,
    }
    return ir


def mock_parameters():
    return {"precision": 8}


def mock_form_ir():
    ir = basic_class_properties("mock_form_classname")
    ir.update({
        "id": 0,
        "prefix": "Foo",
        "signature": "mock form signature",
        "rank": 2,
        "num_coefficients": 3,
        "original_coefficient_position": [0, 2],
    })

    ir.update({
        "create_coordinate_finite_element": ["mock_coordinate_finite_element_classname"],
        "create_coordinate_dofmap": ["mock_coordinate_dofmap_classname"],
        "create_coordinate_mapping": ["mock_coordinate_mapping_classname"],
        "create_finite_element": ["mock_finite_element_classname_%d" % (i,) for i in range(ir["num_coefficients"])],
        "create_dofmap": ["mock_dofmap_classname_%d" % (i,) for i in range(ir["num_coefficients"])],
    })

    # These are the method names in ufc::form that are specialized for each integral type
    template = "max_%s_subdomain_id"
    for i, integral_type in enumerate(ufc_integral_types):
        key = template % integral_type
        ir[key] = i  # just faking some integers

    template = "has_%s_integrals"
    for i, integral_type in enumerate(ufc_integral_types):
        key = template % integral_type
        ir[key] = (i % 2 == 0)  # faking some bools

    template = "create_%s_integral"
    for i, integral_type in enumerate(ufc_integral_types):
        key = template % integral_type
        subdomain_ids = list(range(i))
        classnames = [key.replace("create_", "") + str(j) for j in subdomain_ids]
        ir[key] = (subdomain_ids, classnames)

    template = "create_default_%s_integral"
    for i, integral_type in enumerate(ufc_integral_types):
        key = template % integral_type
        classname = key.replace("create_", "")
        ir[key] = classname

    return ir


def mock_dofmap_ir():
    # Note: The mock data here is not necessarily consistent,
    # it just needs to have the right types to exercise the code generation.
    ir = basic_class_properties("mock_dofmap_classname")
    num_dofs_per_entity = [1, 1, 1]
    entity_dofs = [[(0,), (1,), (2,)],
                   [(3,), (4,), (5,)],
                   [(6,)]]
    entity_closure_dofs = {
        (0, 0): [0],
        (0, 1): [1],
        (0, 2): [2],
        (1, 0): [0, 1],
        (1, 1): [0, 2],
        (1, 2): [1, 2],
        (2, 0): [0, 1, 2],
        }
    ir.update({
        "signature": "mock element signature",
        "geometric_dimension": 3,
        "topological_dimension": 2,
        "global_dimension": ([3, 2, 1], 4),
        "tabulate_dofs": ([[[(0,), (1,), (2,)], [(3,), (4,), (5,)], [(6,)]], None], [7, 1], True, [False, True]),
        "tabulate_facet_dofs": [[0, 1, 2], [1, 2, 3], [0, 2, 3]],
        "tabulate_entity_dofs": (entity_dofs, num_dofs_per_entity),
        "tabulate_entity_closure_dofs": (entity_closure_dofs, entity_dofs, num_dofs_per_entity),
        "needs_mesh_entities": [True, False, True],
        "num_global_support_dofs": 1,
        "num_element_support_dofs": 6,
        "num_element_dofs": 7,
        "num_entity_dofs": num_dofs_per_entity,
        "num_entity_closure_dofs": [3, 6, 10],
        "num_facet_dofs": 7,
        "num_sub_dofmaps": 3,
        "create_sub_dofmap": ["mock_dofmap_classname_sub_%d" % (i,) for i in range(3)],
    })
    return ir


def mock_evaluate_basis_ir():
    dofs_data = [
        {
            "embedded_degree": 5,
            "num_components": 2,
            "num_expansion_members": 7,
            "coeffs": [list(range(20, 27)), list(range(30, 37))],
            "reference_offset": 3,
            "physical_offset": 5,
            "mapping": "affine",
            "dmats": numpy.zeros([]), # FIXME: Mock dmats data
        },
        {
            "embedded_degree": 1,
            "num_components": 2,
            "num_expansion_members": 7,
            "coeffs": [list(range(7)), list(range(10, 17))],
            "reference_offset": 3,
            "physical_offset": 5,
            "mapping": "affine",
            "dmats": numpy.zeros([]), # FIXME: Mock dmats data
        }
    ]
    data = {
        "cellname": "triangle",
        "geometric_dimension": 3,
        "topological_dimension": 2,
        "reference_value_size": 3,
        "physical_value_size": 2,
        "dofs_data": dofs_data,
        "needs_oriented": True,
        "space_dimension": 14,
        "max_degree": 5,
    }
    return data


@pytest.mark.skipif(True, reason="mock ir for evaluate basis is currently incomplete")
def test_mock_evaluate_basis():
    from ffc.uflacs.backends.ufc.evaluatebasis import generate_evaluate_reference_basis
    import ffc.uflacs.language.cnodes as L
    parameters = mock_parameters()
    data = mock_evaluate_basis_ir()  # FIXME: mock data is currently incomplete
    code = generate_evaluate_reference_basis(L, data, parameters)
    print(code)


def mock_finite_element_ir():
    # Note: The mock data here is not necessarily consistent,
    # it just needs to have the right types to exercise the code generation.
    ir = basic_class_properties("mock_finite_element_classname")
    ebir = mock_evaluate_basis_ir()
    ir.update({
        "signature": "mock element signature",
        "cell_shape": "mock_cell_shape",
        "geometric_dimension": 3,
        "topological_dimension": 2,
        "degree": 2,
        "family": "Lagrange",
        "value_dimension": (3, 3),
        "reference_value_dimension": (2, 2),
        "space_dimension": 6,
        "tabulate_dof_coordinates": {"gdim": 3, "tdim": 2, "points": [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)]},
        "evaluate_reference_basis": ebir,
        "evaluate_reference_basis_derivatives": ebir,
        "evaluate_basis": ebir,
        "evaluate_basis_derivatives": ebir,
        "evaluate_basis_all": ebir,
        "evaluate_basis_derivatives_all": ebir,
        "evaluate_dof": "fixme",
        "evaluate_dofs": "fixme",
        "interpolate_vertex_values": "fixme",
        "num_sub_elements": 3,
        "create_sub_element": ["mock_finite_element_classname_sub_%d" % (i,) for i in range(3)],
    })
    return ir


def mock_integral_ir():
    # Note: The mock data here is not necessarily consistent,
    # it just needs to have the right types to exercise the code generation.
    ir = basic_class_properties("mock_integral_classname")
    ir.update({
        "representation": "uflacs",
        "integral_metadata": [{"dummy": "integral metadata"}],
        "integrals_metadata": {"dummy": "integrals metadata"},
        "prefix": "Foo",
        "enabled_coefficients": [True, False, True],
        "tabulate_tensor": "    mock_body_of_tabulate_tensor();",
        "num_cells": 1,
    })
    return ir


def mock_coordinate_mapping_ir():
    ir = basic_class_properties("mock_coordinate_mapping_classname")
    ir.update({
        "signature": "mock_coordinate_mapping_signature",
        "cell_shape": "mock_cell_shape",
        "geometric_dimension": 3,
        "topological_dimension": 2,
        "create_coordinate_finite_element": "mock_coordinate_finite_element_classname",
        "create_coordinate_dofmap": "mock_coordinate_dofmap_classname",
        "scalar_coordinate_finite_element_classname": "mock_scalar_coordinate_finite_element_classname",
        "coordinate_element_degree": 2,
        "num_scalar_coordinate_element_dofs": 5,
        "tables": {"x0": numpy.ones((5,)),
                   "xm": numpy.ones((5,)),
                   "J0": numpy.ones((2, 5)),
                   "Jm": numpy.ones((2, 5)), }
    })
    return ir


def compile_mock_coordinate_mapping():
    parameters = mock_parameters()
    ir = mock_coordinate_mapping_ir()
    gen = ufc_coordinate_mapping()
    return gen.generate(L, ir, parameters)


def compile_mock_form():
    parameters = mock_parameters()
    ir = mock_form_ir()
    gen = ufc_form()
    return gen.generate(L, ir, parameters)


def compile_mock_dofmap():
    parameters = mock_parameters()
    ir = mock_dofmap_ir()
    gen = ufc_dofmap()
    return gen.generate(L, ir, parameters)


def compile_mock_finite_element():
    parameters = mock_parameters()
    ir = mock_finite_element_ir()
    gen = ufc_finite_element()
    return gen.generate(L, ir, parameters)


def compile_mock_integral(integral_type):
    parameters = mock_parameters()
    ir = mock_integral_ir()
    gen = eval("ufc_%s_integral" % integral_type)()
    return gen.generate(L, ir, parameters)


def compile_mock_all():
    mocks = [compile_mock_integral(integral_type) for integral_type in ufc_integral_types]
    mocks += [compile_mock_form(), compile_mock_dofmap(), compile_mock_finite_element()]
    return '\n\n'.join(mocks)


def test_mock_coordinate_mapping():
    h, cpp = compile_mock_coordinate_mapping()
    print(h)
    print(cpp)


def test_mock_form():
    h, cpp = compile_mock_form()
    print(h)
    print(cpp)


def test_mock_dofmap():
    h, cpp = compile_mock_dofmap()
    print(h)
    print(cpp)


@pytest.mark.skipif(True, reason="mock ir for evaluate basis is currently incomplete")
def test_mock_finite_element():
    h, cpp = compile_mock_finite_element()
    print(h)
    print(cpp)


def test_mock_integral():
    for integral_type in ufc_integral_types:
        h, cpp = compile_mock_integral(integral_type)
        print(h)
        print(cpp)


def test_foo_integral_properties():
    parameters = mock_parameters()
    ir = mock_form_ir()
    assert "cell_integral" in ufc_form.create_cell_integral.__doc__
    assert "return" in str(ufc_form().create_cell_integral(L, ir, parameters))


def test_mock_extract_function():
    h, cpp = compile_mock_coordinate_mapping()
    name = "compute_reference_coordinates"
    print("/// Extracted", name, ":")
    print("/// begin")
    print(extract_function(name, cpp))
    print("/// end")


def test_debug_by_printing_extracted_function():
    h, cpp = compile_mock_coordinate_mapping()
    # name = "compute_reference_coordinates"
    # name = "compute_physical_coordinates"
    # name = "compute_jacobians"
    name = "compute_jacobian_inverses"
    print("/// Extracted", name, ":")
    print("/// begin")
    print(extract_function(name, cpp))
    print("/// end")


"""
Missing:

finite_element:
evaluate_basis*
evaluate_dof
interpolate_vertex_values
Ok? tabulate_dof_coordinates

integrals:
tabulate_tensor

all:
everything with classnames
"""
