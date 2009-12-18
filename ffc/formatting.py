"""
This module implements the formatting of UFC code from a given
dictionary of generated C++ code for the body of each UFC function.

It relies on templates for UFC code available as part of the module
ufc_utils. It also relies on templates for generation of specialized
DOLFIN additions available as part of the module dolfin_utils.
"""

__author__ = "Anders Logg (logg@simula.no) and friends"
__date__ = "2009-12-16"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-18

# Python modules
import os

# UFC finite_element templates
from ufc_utils import finite_element_combined
from ufc_utils import finite_element_header
from ufc_utils import finite_element_implementation

# UFC dof_map templates
from ufc_utils import dof_map_combined
from ufc_utils import dof_map_header
from ufc_utils import dof_map_implementation

# FFC modules
from ffc import codesnippets
from ffc.log import info, error
from ffc.constants import FFC_VERSION

def format_ufc(codes, prefix, options):
    "Format given code in UFC format."

    # Generate code for header
    output = _generate_header(prefix, options)

    # Iterate over codes
    for (i, code) in enumerate(codes):

        # Extract generated code
        code_forms, code_elements, code_dofmaps = code

        # Generate code for elements
        for code_element in code_elements:
            output += finite_element_combined % code_element

        # Generate code for dofmaps
        for code_dofmap in code_dofmaps:
            output += dof_map_combined % code_dofmap

    # Generate code for footer
    output += _generate_footer()

    # Write generated code to file
    _write_file(output, prefix, options)

def _write_file(output, prefix, options):
    "Write generated code to file."
    prefix = prefix.split(os.path.join(' ',' ').split()[0])[-1]
    full_prefix = os.path.join(options["output_dir"], prefix)
    filename = "%s.h" % full_prefix
    file = open(filename, "w")
    file.write(output)
    file.close()
    info("Output written to " + filename + ".")

def _generate_header(prefix, options):
    "Generate code for header."
    args = {"version": FFC_VERSION, "prefix_upper": prefix.upper()}
    if options["format"] == "ufc":
        return codesnippets.header_ufc % args
    elif options["format"] == "dolfin":
        return codesnippets.header_dolfin % args
    else:
        error("Unable to format code, unknown format \"%s\".", options["format"])

def _generate_footer():
    "Generate code for footer."
    return codesnippets.footer
