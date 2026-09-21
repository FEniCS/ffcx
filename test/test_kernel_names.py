# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import re

import basix.ufl
import pytest
import ufl

from ffcx.compiler import compile_ufl_objects
from ffcx.options import get_options


def create_mass_form(metadata=None):
    """Create a P1 mass form."""
    coordinate_element = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
    domain = ufl.Mesh(coordinate_element)
    element = basix.ufl.element("Lagrange", "triangle", 1)
    space = ufl.FunctionSpace(domain, element)
    u = ufl.TrialFunction(space)
    v = ufl.TestFunction(space)
    return ufl.inner(u, v) * ufl.dx(metadata=metadata)


def compile_form(form, name="mass", namespace="demo", scalar_type="float64"):
    """Compile a named form and return its header and implementation."""
    options = get_options({"scalar_type": scalar_type})
    code, _ = compile_ufl_objects(
        [form], options=options, object_names={id(form): name}, namespace=namespace
    )
    return code


def test_kernel_name_from_ufl_object_name():
    """Kernel names are derived from form names in UFL files."""
    header, implementation = compile_form(create_mass_form())
    kernel_name = "tabulate_tensor_demo_mass_cell_0_triangle"
    assert f"extern ufcx_tabulate_tensor_float64 {kernel_name};" in header
    assert f"void {kernel_name}(" in implementation
    assert f".tabulate_tensor_float64 = {kernel_name}," in implementation

    # Retain the old hash-based function as a forwarding symbol.
    assert re.search(
        rf"void tabulate_tensor_integral_[0-9a-f]{{40}}_triangle\([^)]*\).*?{kernel_name}\(",
        implementation,
        flags=re.DOTALL,
    )


def test_kernel_name_from_integral_metadata():
    """Integral metadata overrides the form-derived kernel name."""
    form = create_mass_form({"ffcx_kernel_name": "p1_mass"})
    header, implementation = compile_form(form, name="ignored")
    kernel_name = "tabulate_tensor_demo_p1_mass_triangle"
    assert f"extern ufcx_tabulate_tensor_float64 {kernel_name};" in header
    assert f"void {kernel_name}(" in implementation


def test_metadata_name_without_object_name():
    """Explicit metadata names do not require a named UFL object."""
    form = create_mass_form({"ffcx_kernel_name": "p1_mass"})
    options = get_options()
    code, _ = compile_ufl_objects([form], options=options, namespace="demo")
    header, _ = code
    assert "tabulate_tensor_demo_p1_mass_triangle" in header


def test_hash_name_without_object_name():
    """API compilation without object names retains hash-based names."""
    options = get_options()
    code, _ = compile_ufl_objects([create_mass_form()], options=options, namespace="demo")
    header, implementation = code
    assert "extern ufcx_tabulate_tensor_float64" not in header
    assert "tabulate_tensor_demo_mass" not in implementation
    assert re.search(r"void tabulate_tensor_integral_[0-9a-f]{40}_triangle\(", implementation)


@pytest.mark.parametrize("name", ["has-hyphen", "contains space", "λ"])
def test_invalid_kernel_name(name):
    """Kernel names must be valid C identifiers."""
    form = create_mass_form({"ffcx_kernel_name": name})
    with pytest.raises(ValueError, match="Invalid FFCx kernel name"):
        compile_form(form)


def test_non_string_kernel_name():
    """Kernel names must be strings."""
    form = create_mass_form({"ffcx_kernel_name": 3})
    with pytest.raises(ValueError, match="must be a string"):
        compile_form(form)


def test_inconsistent_kernel_name_metadata():
    """Integrals combined into one kernel must use one explicit name."""
    form = create_mass_form()
    integral = form.integrals()[0]
    dx0 = ufl.Measure("dx", domain=integral.ufl_domain(), metadata={"ffcx_kernel_name": "mass_0"})
    dx1 = ufl.Measure("dx", domain=integral.ufl_domain(), metadata={"ffcx_kernel_name": "mass_1"})
    form = integral.integrand() * dx0 + integral.integrand() * dx1
    with pytest.raises(ValueError, match="must use the same 'ffcx_kernel_name'"):
        compile_form(form)


def test_duplicate_kernel_names():
    """Explicit kernel names must be unique within generated code."""
    form0 = create_mass_form({"ffcx_kernel_name": "mass"})
    form1 = create_mass_form({"ffcx_kernel_name": "mass"})
    options = get_options()
    with pytest.raises(ValueError, match="Duplicate FFCx kernel name"):
        compile_ufl_objects(
            [form0, form1],
            options=options,
            object_names={id(form0): "a", id(form1): "b"},
            namespace="demo",
        )
