"""
Compiler stage 5: Code formatting
---------------------------------

This module implements the formatting of UFC code from a given
dictionary of generated C++ code for the body of each UFC function.

It relies on templates for UFC code available as part of the module
ufc_utils.
"""

# Copyright (C) 2009 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2009-12-16
# Last changed: 2011-02-21

# Python modules
import os

# UFC code generation templates
from ufc_utils import templates

# FFC modules
from ffc.log import info, error, begin, end, dstr
from ffc.constants import FFC_VERSION, UFC_VERSION
from ffc.cpp import format

def format_code(code, wrapper_code, prefix, parameters):
    "Format given code in UFC format."

    begin("Compiler stage 5: Formatting code")

    # Extract code
    code_elements, code_dofmaps, code_integrals, code_forms = code

    # Header and implementation code
    code_h = ""
    code_c = ""

    # Generate code for comment on top of file
    code_h += _generate_comment(parameters) + "\n"
    code_c += _generate_comment(parameters) + "\n"

    # Generate code for header
    code_h += format["header_h"] % {"prefix_upper": prefix.upper()}
    code_h += _generate_additional_includes(code_integrals)  + "\n"
    code_c += format["header_c"] % {"prefix": prefix}

    # Generate code for elements
    if code_elements:
        for code_element in code_elements:
            code_h += _format_h("finite_element", code_element, parameters)
            code_c += _format_c("finite_element", code_element, parameters)

    # Generate code for dofmaps
    if code_dofmaps:
        for code_dofmap in code_dofmaps:
            code_h += _format_h("dofmap", code_dofmap, parameters)
            code_c += _format_c("dofmap", code_dofmap, parameters)

    # Generate code for integrals
    if code_integrals:
        for code_integral in code_integrals:
            if "cell_integral" in code_integral["classname"]:
                code_h += _format_h("cell_integral", code_integral, parameters)
                code_c += _format_c("cell_integral", code_integral, parameters)
            elif "exterior_facet_integral" in code_integral["classname"]:
                code_h += _format_h("exterior_facet_integral", code_integral, parameters)
                code_c += _format_c("exterior_facet_integral", code_integral, parameters)
            elif "interior_facet_integral" in code_integral["classname"]:
                code_h += _format_h("interior_facet_integral", code_integral, parameters)
                code_c += _format_c("interior_facet_integral", code_integral, parameters)
            elif "point_integral" in code_integral["classname"]:
                code_h += _format_h("point_integral", code_integral, parameters)
                code_c += _format_c("point_integral", code_integral, parameters)

    # Generate code for form
    if code_forms:
        for code_form in code_forms:
            code_h += _format_h("form", code_form, parameters)
            code_c += _format_c("form", code_form, parameters)

    # Add wrappers
    if wrapper_code:
        code_h += wrapper_code

    # Generate code for footer
    code_h += format["footer"]

    # Write file(s)
    if parameters["split"]:
        _write_file(code_h, prefix, ".h", parameters)
        _write_file(code_c, prefix, ".cpp", parameters)
    else:
        _write_file(code_h, prefix, ".h", parameters)

    end()

def _format_h(class_type, code, parameters):
    "Format header code for given class type."
    if parameters["split"]:
        return templates[class_type + "_header"] % code + "\n"
    else:
        return templates[class_type + "_combined"] % code + "\n"

def _format_c(class_type, code, parameters):
    "Format implementation code for given class type."
    if parameters["split"]:
        return templates[class_type + "_implementation"] % code + "\n"
    else:
        return ""

def _write_file(output, prefix, postfix, parameters):
    "Write generated code to file."
    prefix = prefix.split(os.path.join(' ',' ').split()[0])[-1]
    full_prefix = os.path.join(parameters["output_dir"], prefix)
    filename = "%s%s" % (full_prefix, postfix)
    hfile = open(filename, "w")
    hfile.write(output)
    hfile.close()
    info("Output written to " + filename + ".")

def _generate_comment(parameters):
    "Generate code for comment on top of file."

    # Generate top level comment
    args = {"ffc_version": FFC_VERSION, "ufc_version": UFC_VERSION}
    if parameters["format"] == "ufc":
        comment = format["ufc comment"] % args
    elif parameters["format"] == "dolfin":
        comment = format["dolfin comment"] % args
    else:
        error("Unable to format code, unknown format \"%s\".", parameters["format"])

    # Add parameter information
    comment += format["comment"]("") + "\n"
    comment += format["comment"]("This code was generated with the following parameters:") + "\n"
    comment += format["comment"]("")
    comment += "\n".join([""] + [format["comment"]("  " + l) for l in dstr(parameters).split("\n")][:-1])
    comment += "\n"

    return comment

def _generate_additional_includes(codes):
    s = set()
    for code in codes:
        if "additional_includes_set" in code:
            s.update(code["additional_includes_set"])
    if s:
        return "\n".join(list(s)) + "\n"
    return ""

