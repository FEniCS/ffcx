"Assignment of names to variables"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2008-10-29 -- 2008-10-29"
__copyright__ = "Copyright (C) 2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
from ffc.common.debug import debug

# FFC language modules
from ffc.compiler.language.algebra import Function

def extract_coefficient_names(global_variables):
    "Assign names (attribute .name) to all coefficients found in global_variables"

    debug("Extracting coefficient names...")

    # Assign coefficient names
    coefficient_names = {}
    d = {}
    for name in global_variables:
        variable = global_variables[name]
        if isinstance(variable, Function):
            coefficient_names[variable.n0.index] = str(name)
    coefficient_names = [coefficient_names[i] for i in range(len(coefficient_names))]

    debug("done")

    # Print coefficient names
    debug("Coefficients: " + ", ".join(coefficient_names))

    return coefficient_names
