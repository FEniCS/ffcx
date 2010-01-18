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

# Last changed: 2010-01-18

# Python modules
import os

# UFC code generation templates
from ufc_utils import finite_element_combined
from ufc_utils import dof_map_combined
from ufc_utils import cell_integral_combined
from ufc_utils import exterior_facet_integral_combined
from ufc_utils import interior_facet_integral_combined
from ufc_utils import form_combined

# FFC modules
from ffc import codesnippets
from ffc.log import info, error, begin, end
from ffc.constants import FFC_VERSION, UFC_VERSION

def format_code(codes, wrapper_code, prefix, options):
    "Format given code in UFC format."

    begin("Compiler stage 5: Formatting code")

    # Generate code for header
    output = _generate_header(prefix, options) + "\n\n"

    # Iterate over codes
    for code_elements, code_dofmaps, code_integrals, code_form in codes:

        # Generate code for elements
        if code_elements:
            for code_element in code_elements:
                output += finite_element_combined % code_element + "\n"

        # Generate code for dofmaps
        if code_dofmaps:
            for code_dofmap in code_dofmaps:
                output += dof_map_combined % code_dofmap + "\n"

        # Generate code for integrals
        if code_integrals:
            for code_integral in code_integrals[0]:
                output += cell_integral_combined % code_integral + "\n"
            for code_integral in code_integrals[1]:
                output += exterior_facet_integral_combined % code_integral + "\n"
            for code_integral in code_integrals[2]:
                output += interior_facet_integral_combined % code_integral + "\n"

        # Generate code for form
        if code_form:
            output += form_combined % code_form

    # Add wrapper code
    if wrapper_code:
        output += "\n" + wrapper_code

    # Generate code for footer
    output += _generate_footer()

    # Write generated code to file
    _write_file(output, prefix, options)

    end()

def _write_file(output, prefix, options):
    "Write generated code to file."
    prefix = prefix.split(os.path.join(' ',' ').split()[0])[-1]
    full_prefix = os.path.join(options["output_dir"], prefix)
    filename = "%s.h" % full_prefix
    hfile = open(filename, "w")
    hfile.write(output)
    hfile.close()
    info("Output written to " + filename + ".")

def _generate_header(prefix, options):
    "Generate code for header."
    args = {"ffc_version": FFC_VERSION, "ufc_version": UFC_VERSION, "prefix_upper": prefix.upper()}
    if options["format"] == "ufc":
        return codesnippets.header_ufc % args
    elif options["format"] == "dolfin":
        return codesnippets.header_dolfin % args
    else:
        error("Unable to format code, unknown format \"%s\".", options["format"])

def _generate_footer():
    "Generate code for footer."
    return codesnippets.footer
