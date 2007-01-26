"Code generation for form"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-27 -- 2007-01-27"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

def generate_form(form, format):
    """Generate dictionary of code for the given dof map according to
    the given format."""

    code = {}

    # Generate code for signature
    code["signature"] = form.signature

    # Generate code for rank
    code["rank"] = "%d" % form.rank

    # Generate code for num_coefficients
    code["num_coefficients"] = "%d" % form.num_coefficients
    
    return code
