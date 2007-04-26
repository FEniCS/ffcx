"Form analysis and extraction of form data"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-04-26 -- 2007-04-26"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC language modules
from ffc.compiler.language.reassignment import *

# FFC analysis modules
from checks import *
from formdata import *
from elementdata import *

def analyze(form):
    "Analyze form and extract form data"
    
    # Check validity of form
    check_form(form)

    # Reassign form indices
    reassign_indices(form)

    # Check validity of form again
    check_form(form)

    # Extract form data
    form_data = FormData(form)

    # Print a short summary
    debug("")
    debug("Rank of form: %d" % form_data.rank)
    debug("Coefficients: %d" % form_data.num_coefficients)
    debug("Arguments:    %d" % form_data.num_arguments)
    debug("Terms:        %d" % form_data.num_terms)
    
    return form_data
