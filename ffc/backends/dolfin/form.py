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

from .functionspace import generate_typedefs, function_space_template

def generate_form(form, prefix, classname):
    """Generate dolfin wrapper code associated with a form including code
    for function spaces used in form and typedefs

    @param form:
        A UFCFormNames instance
    @param prefix
        Prefix (namespace) added to names.
    @param classname
        Name of Form class.

    """

    # Generate code for "Form_x_FunctionSpace_y" factories
    wrap = function_space_template
    blocks = [wrap.format(**{"prefix" : prefix,
                     "classname" : "{}_FunctionSpace_{}".format(classname, i),
                     "finite_element_classname" : form.ufc_finite_element_classnames[i],
                     "dofmap_classname" : form.ufc_dofmap_classnames[i],
                     "coordinate_map_classname" : form.ufc_coordinate_mapping_classnames[i]}) for i in range(form.rank)]

    # Generate code for "Form_x_FunctionSpace_y" factories
    # wrap = apply_function_space_template
    # blocks += [wrap("{}_FunctionSpace_{}".format(classname, i),
    #                 form.ufc_finite_element_classnames[i],
    #                 form.ufc_dofmap_classnames[i],
    #                 form.ufc_coordinate_mapping_classnames[i]) for i in range(form.rank)]

    # Add factory function typedefs, e.g. Form_L_FunctionSpace_1_factory = CoefficientSpace_f_factory
    template = "constexpr dolfin_function_space_factory_ptr {0}{1}_FunctionSpace_{2} = {0}CoefficientSpace_{3};"
    blocks += [template.format(prefix, classname, form.rank + i, form.coefficient_names[i])
               for i in range(form.num_coefficients)]

    blocks += ["\n"]

    # Generate Form subclass
    blocks += [generate_form_class(form, prefix, classname)]

    # Return code
    return "\n".join(blocks)


def generate_form_class(form, prefix, classname):
    "Generate dolfin wrapper code for a single Form class."

    # Generate data for coefficient assignments
    (number, name) = generate_coefficient_map_data(form)

    # Generate typedefs for FunctionSpace subclasses for Coefficients
    typedefs = ["// Typedefs (function spaces for {})".format(classname),
                generate_typedefs(form, prefix, classname)]

    # Member variables for coefficients
    # members = ["  dolfin::function::CoefficientAssigner %s;" % coefficient
    #            for coefficient in form.coefficient_names]
    members = []
    #typedefs = []

    # Group typedefs and members together for inserting into template
    additionals = "\n".join(typedefs)

    # Wrap functions in class body
    code = apply_form_template(prefix, classname, form, number, name,
                                additionals)
    return code


def generate_coefficient_map_data(form):
    """Generate data for code for the functions Form::coefficient_number
    and Form::coefficient_name."""

    # Handle case of no coefficients
    if form.num_coefficients == 0:
        num = "  return -1;"
        name = "  return NULL;"
        return (num, name)

    # Otherwise create switch
    ifstr = "if "
    num = ""
    name = '  switch (i)\n  {\n'
    for i, coeff in enumerate(form.coefficient_names):
        num += '  %s(strcmp(name, "%s") == 0)\n    return %d;\n' % (ifstr, coeff, i)
        name += '  case %d:\n    return "%s";\n' % (i, coeff)
        ifstr = 'else if '

    num += "\n  return -1;"
    name += "  }\n  return NULL;"

    return (num, name)


form_class_template = """\
// Return the number of the coefficient with this name
int {prefix}{classname}_coefficient_number(const char* name)
{{
{coefficient_number}
}}

// Return the name of the coefficient with this number
const char* {prefix}{classname}_coefficient_name(int i)
{{
{coefficient_name}
}}

dolfin_form* {prefix}{classname}()
{{
  dolfin_form* form = (dolfin_form*) malloc(sizeof(*form));
  form->form = create_{ufc_form};
  form->coefficient_name_map = {prefix}{classname}_coefficient_name;
  form->coefficient_number_map = {prefix}{classname}_coefficient_number;
  return form;
}}

{typedefs}
"""


def apply_form_template(prefix, classname, form, number, name, typedefs):
    args = {"prefix": prefix,
            "classname": classname,
            "ufc_form": form.ufc_form_classname,
            "coefficient_number": number,
            "coefficient_name": name,
            "typedefs": typedefs}
    return form_class_template.format(**args)
