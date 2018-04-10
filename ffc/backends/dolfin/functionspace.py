# -*- coding: utf-8 -*-
# Copyright (C) 2011 Marie E. Rognes
#
# This file is part of FFV (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Based on original implementation by Martin Alnes and Anders Logg

from .includes import snippets


def extract_coefficient_spaces(forms):
    """Extract a list of tuples

      (classname, finite_element_classname, dofmap_classname, coordinate_mapping_classname)

    for the coefficient spaces in the set of given forms. This can
    then be used for input to the function space template."""

    # Extract data for each coefficient space
    spaces = {}
    for form in forms:
        for (i, name) in enumerate(form.coefficient_names):
            # Skip if already considered
            if name in spaces:
                continue

            # Map element name, dof map name, etc to this coefficient
            spaces[name] = ("CoefficientSpace_{}".format(name),
                            form.ufc_finite_element_classnames[form.rank + i],
                            form.ufc_dofmap_classnames[form.rank + i],
                            form.ufc_coordinate_mapping_classnames[form.rank + i])

    # Return coefficient spaces sorted alphabetically by coefficient name
    names = sorted(spaces.keys())
    return [spaces[name] for name in names]


def generate_typedefs(form, classname):
    """Generate typedefs for test, trial and coefficient spaces relative
    to a function space.

    """

    # Add convenience pointers to factory functions
    template0 = "  constexpr dolfin_function_space_factory_ptr {}_factory = {}_FunctionSpace_{}_factory;"
    factory0 = "\n".join(template0.format(
        snippets["functionspace"][i], classname, i) for i in range(form.rank))

    template1 = "  constexpr dolfin_function_space_factory_ptr CoefficientSpace_{}_factory = {}_FunctionSpace_{}_factory;"
    factory1 = "\n".join(template1.format(
        form.coefficient_names[i], classname, form.rank + i) for i in range(form.num_coefficients))

    code = "\n" + factory0 + "\n" + factory1
    return code


function_space_template = """\
dolfin_function_space* %(classname)s_factory()
{
  /* dolfin_function_space* space = malloc(sizeof(*space)); // In C rather than C++: */
  dolfin_function_space* space = (dolfin_function_space*) malloc(sizeof(*space));
  space->element = create_%(ufc_finite_element_classname)s;
  space->dofmap = create_%(ufc_dofmap_classname)s;
  space->coordinate_mapping = create_%(ufc_coordinate_mapping_classname)s;
  return space;
}

"""


def apply_function_space_template(name, element_name, dofmap_name,
                                  coordinate_mapping):
    args = {"classname": name,
            "ufc_finite_element_classname": element_name,
            "ufc_dofmap_classname": dofmap_name,
            "ufc_coordinate_mapping_classname": coordinate_mapping}
    return function_space_template % args
