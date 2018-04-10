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

from .functionspace import (extract_coefficient_spaces, apply_function_space_template)
from .form import generate_form
from .capsules import UFCElementNames

# NB: generate_dolfin_namespace(...) assumes that if a coefficient has
# the same name in multiple forms, it is indeed the same coefficient:
parameters = {"use_common_coefficient_names": True}


def generate_dolfin_code(prefix, forms, common_function_space=False):
    """Generate complete dolfin wrapper code with given generated names.

    @param prefix:
        String, prefix for all form names
    @param forms:
        List of UFCFormNames instances or single UFCElementNames.
    @param common_function_space:
        True if common function space, otherwise False
    """

    # Tag
    code = "// DOLFIN helper functions wrappers\n"

    # Typedefs for convenience factory functions
    factory_typedefs = """
// Typedefs for convenience factory functions
typedef dolfin_function_space* (*dolfin_function_space_factory_ptr)(void);
typedef dolfin_form* (*dolfin_form_factory_ptr)(void);

"""
    code += factory_typedefs

    # Generate body of dolfin wappers
    if isinstance(forms, UFCElementNames):
        namespace = generate_single_function_space(prefix, forms)
    else:
        namespace = generate_dolfin_namespace(prefix, forms, common_function_space)

    code += "\n".join([namespace])

    return code


def generate_dolfin_namespace(prefix, forms, common_function_space=False):

    # Extract (common) coefficient spaces
    assert(parameters["use_common_coefficient_names"])
    spaces = extract_coefficient_spaces(forms)

    # Generate code for common coefficient spaces
    code = [apply_function_space_template(*space) for space in spaces]

    # Generate code for forms, including function spaces for test/trial functions
    code += [generate_form(form, "Form_{}".format(form.name)) for form in forms]

    # Generate 'top-level' typedefs (Bilinear/Linear & Test/Trial/Function)
    code += [generate_namespace_typedefs(forms, common_function_space)]

    # Wrap code in namespace block
    code = "\n".join(code)

    # Return code
    return code


def generate_single_function_space(prefix, space):
    code = apply_function_space_template("FunctionSpace",
                                         space.ufc_finite_element_classnames[0],
                                         space.ufc_dofmap_classnames[0],
                                         space.ufc_coordinate_mapping_classnames[0])
    return code


def generate_namespace_typedefs(forms, common_function_space):

    # Generate typedefs as (fro, to) pairs of strings
    pairs = []

    # Add typedef for Functional/LinearForm/BilinearForm if only one
    # is present of each
    aliases = ["Functional", "LinearForm", "BilinearForm"]
    extra_aliases = {"LinearForm": "ResidualForm",
                     "BilinearForm": "JacobianForm"}
    for rank in sorted(range(len(aliases)), reverse=True):
        forms_of_rank = [form for form in forms if form.rank == rank]
        if len(forms_of_rank) == 1:
            pairs += [("Form_{}".format(forms_of_rank[0].name), aliases[rank])]
            if aliases[rank] in extra_aliases:
                extra_alias = extra_aliases[aliases[rank]]
                pairs += [("Form_{}".format(forms_of_rank[0].name), extra_alias)]

    # Combine data to typedef code
    typedefs = "\n".join("constexpr dolfin_form_factory_ptr {}_factory = {}_factory;".format(fro, to) for (to, fro) in pairs)

    # Keepin' it simple: Add typedef for function space factory if term applies
    if common_function_space:
        for i, form in enumerate(forms):
            if form.rank:
                # FIXME: Is this naming robust?
                typedefs += "\n\nstatic constexpr dolfin_function_space_factory_ptr FunctionSpace_factory = Form_{}_FunctionSpace_0_factory;".format(form.name)
                break

    # Return typedefs or ""
    if not typedefs:
        return ""
    return "// High-level typedefs\n" + typedefs + "\n"
