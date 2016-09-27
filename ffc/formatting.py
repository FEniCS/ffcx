"""
Compiler stage 5: Code formatting
---------------------------------

This module implements the formatting of UFC code from a given
dictionary of generated C++ code for the body of each UFC function.

It relies on templates for UFC code available as part of the module
ufc_utils.
"""

# Copyright (C) 2009-2015 Anders Logg
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

# Python modules
import os

# FFC modules
from ffc.log import info, error, begin, end, dstr
from ffc import __version__ as FFC_VERSION
from ffc.backends.ufc import __version__ as UFC_VERSION
from ffc.cpp import format, make_classname
from ffc.backends.ufc import templates, visibility_snippet, factory_decl, factory_impl
from ffc.parameters import compilation_relevant_parameters


def format_code(code, wrapper_code, prefix, parameters, jit=False):
    "Format given code in UFC format. Returns two strings with header and source file contents."

    begin("Compiler stage 5: Formatting code")

    # Extract code
    code_elements, code_dofmaps, code_coordinate_mappings, code_integrals, code_forms = code

    # Generate code for comment on top of file
    code_h_pre = _generate_comment(parameters) + "\n"
    code_c_pre = _generate_comment(parameters) + "\n"

    # Generate code for header
    code_h_pre += format["header_h"] % {"prefix_upper": prefix.upper()}
    code_c_pre += format["header_c"] % {"prefix": prefix}

    # Add includes
    includes_h, includes_c = _generate_additional_includes(code)
    code_h_pre += includes_h
    code_c_pre += includes_c

    # Header and implementation code
    code_h = ""
    code_c = ""

    if jit:
        code_c += visibility_snippet

    # Generate code for elements
    for code_element in code_elements:
        code_h += _format_h("finite_element", code_element, parameters, jit)
        code_c += _format_c("finite_element", code_element, parameters, jit)

    # Generate code for dofmaps
    for code_dofmap in code_dofmaps:
        code_h += _format_h("dofmap", code_dofmap, parameters, jit)
        code_c += _format_c("dofmap", code_dofmap, parameters, jit)

    # Generate code for coordinate_mappings
    code_coordinate_mappings = []  # FIXME: This disables output of generated coordinate_mapping class, until implemented properly
    for code_coordinate_mapping in code_coordinate_mappings:
        code_h += _format_h("coordinate_mapping", code_coordinate_mapping, parameters)
        code_c += _format_c("coordinate_mapping", code_coordinate_mapping, parameters)

    # Generate code for integrals
    for code_integral in code_integrals:
        code_h += _format_h(code_integral["class_type"], code_integral, parameters, jit)
        code_c += _format_c(code_integral["class_type"], code_integral, parameters, jit)

    # Generate code for form
    for code_form in code_forms:
        code_h += _format_h("form", code_form, parameters, jit)
        code_c += _format_c("form", code_form, parameters, jit)

    # Add factory functions named "..._main" to construct the
    # main element and dofmap in this module
    if jit and not code_forms:
        kind = "finite_element"
        pub = make_classname(prefix, kind, "main")
        code_h += factory_decl % {
            "basename": "ufc::%s" % kind,
            "publicname": pub,
        }
        code_c += factory_impl % {
            "basename": "ufc::%s" % kind,
            "publicname": pub,
            "privatename": code_elements[-1]["classname"]
        }
        kind = "dofmap"
        pub = make_classname(prefix, kind, "main")
        code_h += factory_decl % {
            "basename": "ufc::%s" % kind,
            "publicname": pub,
        }
        code_c += factory_impl % {
            "basename": "ufc::%s" % kind,
            "publicname": pub,
            "privatename": code_dofmaps[-1]["classname"]
        }

    # Add wrappers
    if wrapper_code:
        code_h += wrapper_code

    # Generate code for footer
    code_h += format["footer"]

    # Add headers to body
    code_h = code_h_pre + code_h
    if code_c:
        code_c = code_c_pre + code_c

    end()

    return code_h, code_c


def write_code(code_h, code_c, prefix, parameters):
    # Write file(s)
    _write_file(code_h, prefix, ".h", parameters)
    if code_c:
        _write_file(code_c, prefix, ".cpp", parameters)


def _format_h(class_type, code, parameters, jit=False):
    "Format header code for given class type."
    if jit:
        return templates[class_type + "_jit_header"] % code + "\n"
    elif parameters["split"]:
        return templates[class_type + "_header"] % code + "\n"
    else:
        return templates[class_type + "_combined"] % code + "\n"


def _format_c(class_type, code, parameters, jit=False):
    "Format implementation code for given class type."
    if jit:
        return templates[class_type + "_jit_implementation"] % code + "\n"
    elif parameters["split"]:
        return templates[class_type + "_implementation"] % code + "\n"
    else:
        return ""


def _write_file(output, prefix, postfix, parameters):
    "Write generated code to file."
    filename = os.path.join(parameters["output_dir"], prefix + postfix)
    with open(filename, "w") as hfile:
        hfile.write(output)
    info("Output written to " + filename + ".")


def _generate_comment(parameters):
    "Generate code for comment on top of file."

    # Drop irrelevant parameters
    parameters = compilation_relevant_parameters(parameters)

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


def _generate_additional_includes(code):
    s = set()
    s.add("#include <ufc.h>")

    for code_foo in code:
        # FIXME: Avoid adding these includes if we don't need them
        # s.add("#include <cmath>")
        s.add("#include <stdexcept>")

        for c in code_foo:
            s.update(c.get("additional_includes_set", ()))

    if s:
        includes = "\n".join(sorted(s)) + "\n"
    else:
        includes = ""

    # TODO: move includes to cpp file if possible
    return includes, ""
