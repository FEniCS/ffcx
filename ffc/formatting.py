"""
Compiler stage 5: Code formatting
---------------------------------

This module implements the formatting of UFC code from a given
dictionary of generated C++ code for the body of each UFC function.

It relies on templates for UFC code available as part of the module
ufc_utils.
"""

__author__ = "Anders Logg (logg@simula.no) and friends"
__date__ = "2009-12-16"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-31

# Python modules
import os

# UFC code generation templates
from ufc_utils import templates

# FFC modules
from ffc.log import info, error, begin, end
from ffc.constants import FFC_VERSION, UFC_VERSION
from ffc.cpp import format

def format_code(code, wrapper_code, prefix, options):
    "Format given code in UFC format."

    begin("Compiler stage 5: Formatting code")

    # Extract code
    code_elements, code_dofmaps, code_integrals, code_forms = code

    # Header and implementation code
    code_h = ""
    code_c = ""

    # Generate code for comment on top of file
    code_h += _generate_comment(options) + "\n"
    code_c += _generate_comment(options) + "\n"

    # Generate code for header
    code_h += format["header"] % {"prefix_upper": prefix.upper()} + "\n"

    # Generate code for elements
    if code_elements:
        for code_element in code_elements:
            code_h += _format_h("finite_element", code_element, options)
            code_c += _format_c("finite_element", code_element, options)

    # Generate code for dofmaps
    if code_dofmaps:
        for code_dofmap in code_dofmaps:
            code_h += _format_h("dof_map", code_dofmap, options)
            code_c += _format_c("dof_map", code_dofmap, options)

    # Generate code for integrals
    if code_integrals:
        for code_integral in code_integrals:
            if "cell_integral" in code_integral["classname"]:
                code_h += _format_h("cell_integral", code_integral, options)
                code_c += _format_c("cell_integral", code_integral, options)
            elif "exterior_facet_integral" in code_integral["classname"]:
                code_h += _format_h("exterior_facet_integral", code_integral, options)
                code_c += _format_c("exterior_facet_integral", code_integral, options)
            elif "interior_facet_integral" in code_integral["classname"]:
                code_h += _format_h("interior_facet_integral", code_integral, options)
                code_c += _format_c("interior_facet_integral", code_integral, options)

    # Generate code for form
    if code_forms:
        for code_form in code_forms:
            code_h += _format_h("form", code_form, options)
            code_c += _format_c("form", code_form, options)

    # Add wrappers
    if wrapper_code:
        code_h += wrapper_code

    # Generate code for footer
    code_h += format["footer"]

    # Write file(s)
    if options["split"]:
        _write_file(code_h, prefix, ".h", options)
        _write_file(code_c, prefix, ".cpp", options)
    else:
        _write_file(code_h, prefix, ".h", options)

    end()

def _format_h(class_type, code, options):
    "Format header code for given class type."
    if options["split"]:
        return templates[class_type + "_header"] % code + "\n"
    else:
        return templates[class_type + "_combined"] % code + "\n"

def _format_c(class_type, code, options):
    "Format implementation code for given class type."
    if options["split"]:
        return templates[class_type + "_implementation"] % code + "\n"
    else:
        return ""

def _write_file(output, prefix, postfix, options):
    "Write generated code to file."
    prefix = prefix.split(os.path.join(' ',' ').split()[0])[-1]
    full_prefix = os.path.join(options["output_dir"], prefix)
    filename = "%s%s" % (full_prefix, postfix)
    hfile = open(filename, "w")
    hfile.write(output)
    hfile.close()
    info("Output written to " + filename + ".")

def _generate_comment(options):
    "Generate code for comment on top of file."
    args = {"ffc_version": FFC_VERSION, "ufc_version": UFC_VERSION}
    if options["format"] == "ufc":
        return format["ufc comment"] % args
    elif options["format"] == "dolfin":
        return format["dolfin comment"] % args
    else:
        error("Unable to format code, unknown format \"%s\".", options["format"])
