"Adding appropriate ranges to undetermined indices."

__author__ = "Marie Rognes (meg@math.uio.no)"
__date__ = "2007-04-30 -- 2007-04-30"
__copyright__ = "Copyright (C) 2007 "
__license__  = "GNU GPL Version 2"

# Modified by Anders Logg 2007

# Python modules
import sys

# FFC common modules
from ffc.common.debug import *
from ffc.common.exceptions import *

# FFC compiler.language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.algebra import *
from ffc.compiler.language.tokens import *

def determine_index_ranges(form):
    "Go through the form and set appropriate index ranges."
    debug("Determining the range of the indices...")
    for monomial in form.monomials:
        # Skip constants for now (will be removed and replaced by DG(0) functions)
        #for constant in monomial.constants:
        #    set_range(constant.number, 0)
        for coefficient in monomial.coefficients:
            # The index corresponds to the summation space dimension
            set_range(coefficient.n0, coefficient.e0.space_dimension())
            set_range(coefficient.n1, coefficient.e1.space_dimension())
            set_range(coefficient.index, coefficient.e1.space_dimension())
        for transform in monomial.transforms:
            # The indices of the transforms correspond to the cell dimension
            set_range(transform.index0, transform.element.cell_dimension())
            set_range(transform.index1, transform.element.cell_dimension())
        for basisfunction in monomial.basisfunctions:
            # The index of the basisfunction corresponds to the space dimension
            set_range(basisfunction.index, basisfunction.element.space_dimension())
            for i in range(len(basisfunction.component)):
                # The component corresponds to the value dimension
                set_range(basisfunction.component[i], basisfunction.element.value_dimension(i))
            for d in basisfunction.derivatives:
                # The index of the derivatives corresponds to the cell dimension
                set_range(d.index, d.element.cell_dimension())
    debug("done")

def set_range(index, dimension):
    "Set index range for indices that don't already have their indices specified"
    if isinstance(index, Index):
        if not index.range: # i.e. index.range == none:
            index.range = range(dimension)
    else:
        raise FormError("Cannot set range for non-indices!")
