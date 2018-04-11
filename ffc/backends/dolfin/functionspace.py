# -*- coding: utf-8 -*-
# Copyright (C) 2011 Marie E. Rognes
#
# This file is part of FFV (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Based on original implementation by Martin Alnes and Anders Logg


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


def generate_typedefs(form, prefix, classname):
    """Generate typedefs for test, trial and coefficient spaces relative
    to a function space.

    """

    snippets = {"functionspace" : ("TestSpace", "TrialSpace")}

    # Add convenience pointers to factory functions
    template0 = "constexpr dolfin_function_space_factory_ptr {0}{2}{1} = {0}{2}_FunctionSpace_{3};"
    factory0 = "\n".join(template0.format(prefix,
        snippets["functionspace"][i], classname, i) for i in range(form.rank))

    template1 = "constexpr dolfin_function_space_factory_ptr {0}{2}CoefficientSpace_{1} = {0}{2}_FunctionSpace_{3};"
    factory1 = "\n".join(template1.format(prefix,
        form.coefficient_names[i], classname, form.rank + i) for i in range(form.num_coefficients))

    code = factory0 + "\n" + factory1
    return code


function_space_template = """\
dolfin_function_space* {prefix}{classname}()
{{
  /* dolfin_function_space* space = malloc(sizeof(*space)); // In C rather than C++: */
  dolfin_function_space* space = (dolfin_function_space*) malloc(sizeof(*space));
  space->element = create_{finite_element_classname};
  space->dofmap = create_{dofmap_classname};
  space->coordinate_mapping = create_{coordinate_map_classname};
  return space;
}}
"""
