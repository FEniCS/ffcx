"Code generation for form"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-27 -- 2007-02-06"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

def generate_form(form_data, format):
    """Generate dictionary of code for the given form data map
    according to the given format"""

    code = {}

    # Generate code for signature
    code["signature"] = form_data.signature

    # Generate code for rank
    code["rank"] = "%d" % form_data.rank

    # Generate code for num_coefficients
    code["num_coefficients"] = "%d" % form_data.num_coefficients
    
    return code
