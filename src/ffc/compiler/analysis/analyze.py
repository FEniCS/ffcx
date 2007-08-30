"Form analysis and extraction of form data"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-04-26 -- 2007-04-26"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
from ffc.common.debug import *

# FFC language modules
from ffc.compiler.language.reassignment import *

# FFC analysis modules
from checks import *
from formdata import *
from elementdata import *
from simplify import *
from ranges import *

def analyze(form):
    "Analyze form and extract form data"

    raw_form = str(form)

    # Check validity of form
    check_form(form)

    # Reassign form indices
    reassign_indices(form)
    reassigned_form = str(form)

    # Check validity of form again
    check_form(form)

    # Simplify form
    simplify(form)
    simplified_form = str(form)

    # Check validity of form again
    check_form(form)

    # Determine range of indices:
    determine_index_ranges(form)

    check_form(form)

    # Extract form data
    form_data = FormData(form)

    # Print form before and after reassignment and simplification
    debug("\nRaw input form:\n", 1)
    debug("    " + raw_form, 1)
    debug("\nAfter index reassignment:\n", 1)
    debug("    " + reassigned_form, 1)
    debug("\nAfter simplification:\n", 1)
    debug("    " + simplified_form, 1)
    
    # Print a short summary
    debug("")
    debug("Rank of form: %d" % form_data.rank)
    debug("Coefficients: %d" % form_data.num_coefficients)
    debug("Arguments:    %d" % form_data.num_arguments)
    debug("Terms:        %d" % form_data.num_terms)
    
    return form_data
