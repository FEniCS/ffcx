# -*- coding: utf-8 -*-
# Copyright (C) 2011 Marie E. Rognes
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
#
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Based on original implementation by Martin Alnes and Anders Logg
#
# Modified by Anders Logg 2015
#
# Last changed: 2016-03-02

from .includes import snippets

__all__ = ["apply_function_space_template",
           "extract_coefficient_spaces",
           "generate_typedefs"]


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

            # Map element name, dof map name etc to this coefficient
            spaces[name] = ("CoefficientSpace_%s" % name,
                            form.ufc_finite_element_classnames[form.rank + i],
                            form.ufc_dofmap_classnames[form.rank + i],
                            form.ufc_coordinate_mapping_classnames[form.rank + i])

    # Return coefficient spaces sorted alphabetically by coefficient
    # name
    names = sorted(spaces.keys())
    return [spaces[name] for name in names]


def generate_typedefs(form, classname):
    """Generate typedefs for test, trial and coefficient spaces relative
    to a function space.

    """

    pairs = []

    # # Generate typedef data for test/trial spaces
    # pairs += [("%s_FunctionSpace_%d" % (classname, i),
    #           snippets["functionspace"][i]) for i in range(form.rank)]

    # # Generate typedefs for coefficient spaces
    # pairs += [("%s_FunctionSpace_%d" % (classname, form.rank + i),
    #            "CoefficientSpace_%s" % form.coefficient_names[i])
    #           for i in range(form.num_coefficients)]

    # Combine data to typedef code
    code = "\n".join("  typedef {} {};".format(to, fro) for (to, fro) in pairs)

    # Add convenience pointers to factory functions
    template0 = "  static constexpr dolfin_function_space_factory_ptr {}_factory = {}_FunctionSpace_{}_factory;"
    factory0 = "\n".join(template0.format(
        snippets["functionspace"][i], classname, i) for i in range(form.rank))

    template1 = "  static constexpr dolfin_function_space_factory_ptr CoefficientSpace_{}_factory = {}_FunctionSpace_{}_factory;"
    factory1 = "\n".join(template1.format(
        form.coefficient_names[i], classname, form.rank + i) for i in range(form.num_coefficients))

    code += "\n" + factory0 + "\n" + factory1
    return code


function_space_template = """\
dolfin_function_space* %(classname)s_factory()
{
  /*
  In C rather than C++:
  dolfin_function_space* space = malloc(sizeof(*space));
  */
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
