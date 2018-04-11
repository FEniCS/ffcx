# -*- coding: utf-8 -*-
# Copyright (C) 2011 Marie E. Rognes
#
# This file is part of FFV (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Based on original implementation by Martin Alnes and Anders Logg

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
        # FIXME: update
        #code_blocks = generate_single_function_space(prefix, forms)
        code = apply_function_space_template("FunctionSpace",
                                         space.ufc_finite_element_classnames[0],
                                         space.ufc_dofmap_classnames[0],
                                         space.ufc_coordinate_mapping_classnames[0])
    else:
        # Generate wrappers for a form
        code_blocks = []

        # FIXME: Convert to dict
        # Extract (common) coefficient spaces
        assert(parameters["use_common_coefficient_names"])
        spaces = extract_coefficient_spaces(forms)

        # Generate dolfin_function_space code for common coefficient spaces
        for space in spaces:
            args = {"prefix" : prefix,
                    "classname" : space[0],
                    "finite_element_classname" : space[1],
                    "dofmap_classname" : space[2],
                    "coordinate_map_classname" : space[3]}
            code_blocks += [FUNCTION_SPACE_TEMPLATE.format(**args)]
        #code = [apply_function_space_template(*space) for space in spaces]

        # Generate code for forms (including function spaces for test/trial functions)
        code_blocks += [generate_form(form, prefix, "Form_{}".format(form.name)) for form in forms]

        # Generate 'top-level' typedefs (Bilinear/Linear & Test/Trial/Function)
        code_blocks += [generate_namespace_typedefs(forms, prefix, common_function_space)]

    code += "\n".join(code_blocks)

    return code


def generate_namespace_typedefs(forms, prefix, common_function_space):

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
    typedefs = "\n".join("constexpr dolfin_form_factory_ptr {0}{1} = {0}{2};".format(prefix, fro, to) for (to, fro) in pairs)

    # Keepin' it simple: Add typedef for function space factory if term applies
    if common_function_space:
        for i, form in enumerate(forms):
            if form.rank:
                # FIXME: Is this naming robust?
                typedefs += "\n\nstatic constexpr dolfin_function_space_factory_ptr {0}FunctionSpace = {0}Form_{1}_FunctionSpace_0;".format(prefix, form.name)
                break

    # Return typedefs or ""
    if not typedefs:
        return ""
    return "// High-level typedefs\n" + typedefs + "\n"


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
    assert(len(form.ufc_coordinate_mapping_classnames) == 1)
    wrap = FUNCTION_SPACE_TEMPLATE
    blocks = [wrap.format(**{"prefix" : prefix,
                             "classname" : "{}_FunctionSpace_{}".format(classname, i),
                             "finite_element_classname" : form.ufc_finite_element_classnames[i],
                             "dofmap_classname" : form.ufc_dofmap_classnames[i],
                             "coordinate_map_classname" : form.ufc_coordinate_mapping_classnames[0]}) for i in range(form.rank)]

    # Add factory function typedefs, e.g. Form_L_FunctionSpace_1_factory = CoefficientSpace_f_factory
    blocks += ["// Coefficient function spaces for form \"{}\"".format(classname)]
    template = "constexpr dolfin_function_space_factory_ptr {0}{1}_FunctionSpace_{2} = {0}CoefficientSpace_{3};"
    blocks += [template.format(prefix, classname, form.rank + i, form.coefficient_names[i])
               for i in range(form.num_coefficients)]
    if form.num_coefficients == 0:
        blocks += ["// None"]
    blocks += [""]

    # Generate Form subclass
    blocks += [generate_form_class(form, prefix, classname)]

    # Return code
    return "\n".join(blocks)


def generate_form_class(form, prefix, classname):
    "Generate dolfin wrapper code for a single Form class."

    # Generate data for coefficient assignments
    (number, name) = generate_coefficient_map_data(form)

    # Generate typedefs for FunctionSpace subclasses for Coefficients
    typedefs = "// Typedefs (function spaces for {})\n".format(classname) + generate_function_space_typedefs(form, prefix, classname)

    # Wrap functions in class body
    args = {"prefix": prefix,
            "classname": classname,
            "ufc_form": form.ufc_form_classname,
            "coefficient_number": number,
            "coefficient_name": name,
            "typedefs": typedefs}

    code = FORM_CLASS_TEMPLATE.format(**args)

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
            assert len(form.ufc_coordinate_mapping_classnames) == 1
            spaces[name] = ("CoefficientSpace_{}".format(name),
                            form.ufc_finite_element_classnames[form.rank + i],
                            form.ufc_dofmap_classnames[form.rank + i],
                            form.ufc_coordinate_mapping_classnames[0])

    # Return coefficient spaces sorted alphabetically by coefficient name
    names = sorted(spaces.keys())
    return [spaces[name] for name in names]


def generate_function_space_typedefs(form, prefix, classname):
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


FORM_CLASS_TEMPLATE = """\
// Return the number of the coefficient with this name. Returns -1 if name does not exist.
int {prefix}{classname}_coefficient_number(const char* name)
{{
{coefficient_number}
}}

// Return the name of the coefficient with this number. Returns NULL if index is out-of-range.
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


FUNCTION_SPACE_TEMPLATE = """\
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
