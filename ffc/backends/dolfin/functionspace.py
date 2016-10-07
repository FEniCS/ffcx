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
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Based on original implementation by Martin Alnes and Anders Logg
#
# Modified by Anders Logg 2015
#
# Last changed: 2016-03-02

from .includes import snippets

__all__ = ["apply_function_space_template",
           "apply_multimesh_function_space_template",
           "extract_coefficient_spaces",
           "generate_typedefs"]


def extract_coefficient_spaces(forms):
    """Extract a list of tuples

      (classname, finite_element_classname, dofmap_classname)

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
                            form.ufc_dofmap_classnames[form.rank + i])

    # Return coefficient spaces sorted alphabetically by coefficient
    # name
    names = sorted(spaces.keys())
    return [spaces[name] for name in names]


def generate_typedefs(form, classname, error_control):
    """Generate typedefs for test, trial and coefficient spaces
    relative to a function space."""

    pairs = []

    # Generate typedef data for test/trial spaces
    pairs += [("%s_FunctionSpace_%d" % (classname, i),
              snippets["functionspace"][i]) for i in range(form.rank)]
    if not error_control:  # FIXME: Issue #91
        pairs += [("%s_MultiMeshFunctionSpace_%d" % (classname, i),
                   snippets["multimeshfunctionspace"][i]) for i in range(form.rank)]

    # Generate typedefs for coefficient spaces
    pairs += [("%s_FunctionSpace_%d" % (classname, form.rank + i),
               "CoefficientSpace_%s" % form.coefficient_names[i])
              for i in range(form.num_coefficients)]

    # Combine data to typedef code
    code = "\n".join("  typedef %s %s;" % (to, fro) for (to, fro) in pairs)
    return code


function_space_template = """\
class %(classname)s: public dolfin::FunctionSpace
{
public:

  // Constructor for standard function space
  %(classname)s(std::shared_ptr<const dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          std::make_shared<const dolfin::FiniteElement>(std::make_shared<%(ufc_finite_element_classname)s>()),
                          std::make_shared<const dolfin::DofMap>(std::make_shared<%(ufc_dofmap_classname)s>(), *mesh))
  {
    // Do nothing
  }

  // Constructor for constrained function space
  %(classname)s(std::shared_ptr<const dolfin::Mesh> mesh, std::shared_ptr<const dolfin::SubDomain> constrained_domain):
    dolfin::FunctionSpace(mesh,
                          std::make_shared<const dolfin::FiniteElement>(std::make_shared<%(ufc_finite_element_classname)s>()),
                          std::make_shared<const dolfin::DofMap>(std::make_shared<%(ufc_dofmap_classname)s>(), *mesh, constrained_domain))
  {
    // Do nothing
  }

};
"""


multimesh_function_space_template = """\
class %(classname)s: public dolfin::MultiMeshFunctionSpace
{
public:

  // Constructor for multimesh function space
  %(classname)s(std::shared_ptr<const dolfin::MultiMesh> multimesh): dolfin::MultiMeshFunctionSpace(multimesh)
  {
    // Create and add standard function spaces
    for (std::size_t part = 0; part < multimesh->num_parts(); part++)
    {
      std::shared_ptr<const dolfin::FunctionSpace> V(new %(single_name)s(multimesh->part(part)));
      add(V);
    }

    // Build multimesh function space
    build();
  }

};
"""


def apply_function_space_template(name, element_name, dofmap_name):
    args = {"classname": name,
            "ufc_finite_element_classname": element_name,
            "ufc_dofmap_classname": dofmap_name}
    return function_space_template % args


def apply_multimesh_function_space_template(name, single_name, element_name,
                                            dofmap_name):
    args = {"classname": name,
            "single_name": single_name,
            "ufc_finite_element_classname": element_name,
            "ufc_dofmap_classname": dofmap_name}
    return multimesh_function_space_template % args
