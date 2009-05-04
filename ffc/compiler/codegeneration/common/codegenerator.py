__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-06 -- 2009-05-04"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.

# FFC common modules
from ffc.common.log import debug

# Code generation modules
from finiteelement import generate_finite_elements
from dofmap import generate_dof_maps
from form import generate_form

def generate_common_code(form_data, format):
    "Generate common form code according to given format."

    code = {}

    # Generate code for finite elements
    code["finite_elements"] = generate_finite_elements(form_data, format)

    # Generate code for dof maps
    code["dof_maps"] = generate_dof_maps(form_data, format)

    # Generate code for form
    debug("Generating code for form...")
    code["form"] = generate_form(form_data, format)
    debug("done")

    return code
