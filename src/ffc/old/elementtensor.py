__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-06 -- 2007-01-23"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# FFC common modules
from debug import *
from constants import *
from exceptions import *

# FFC compiler.language modules
from integral import *

# FFC compiler modules
from tensorrepresentation import *


class ElementTensor:
    """An ElementTensor represents the element tensor of a
    multi-linear form and consists of a list of Terms, each containing
    a pair of a ReferenceTensor and a GeometryTensor.

    Attributes:

        terms   - a list of Terms (products A0 * GK)
        aK      - a list of precomputed element tensor declarations
        a0      - a list of precomputed reference tensor declarations
        gK      - a list of precomputed geometry tensor declarations
        num_ops - number of operations in computation of element tensor
    """

    def __init__(self, sum, format, cK_used, gK_used, options):
        "Create ElementTensor."

        # Reset number of operations
        self.num_ops = 0

        # Compute terms
        self.terms = compute_terms(sum, Integral.CELL, None, None, None)
        if len(self.terms) == 0:
            return

        # Compute element tensor declarations
        (self.aK, self.num_ops) = compute_element_tensor(self.terms, format, options)

        # Compute reference tensor declarations
        self.a0 = compute_reference_tensor(self.terms, format)
        
        # Compute geometry tensor declarations
        check_used(self.terms, format, gK_used)
        self.gK = compute_geometry_tensor(self.terms, format, gK_used, cK_used)
