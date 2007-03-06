"Code generation for form"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-27 -- 2007-03-06"
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

    # Generate code for num_cell_integrals
    code["num_cell_integrals"] = "%d" % form_data.num_cell_integrals

    # Generate code for num_exterior_facet_integrals
    code["num_exterior_facet_integrals"] = "%d" % form_data.num_exterior_facet_integrals
    
    # Generate code for num_interior_facet_integrals
    code["num_interior_facet_integrals"] = "%d" % form_data.num_interior_facet_integrals
    
    return code
