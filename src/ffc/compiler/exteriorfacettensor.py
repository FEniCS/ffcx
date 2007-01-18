__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2007-01-180"
__copyright__ = "Copyright (C) 2006-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *
from ffc.common.exceptions import *

# FFC formlanguage modules
from ffc.formlanguage.integral import *

# FFC compiler modules
from tensorrepresentation import *

class ExteriorFacetTensor:
    """An ExteriorElementTensor represents the exterior facet tensor for
    a multi-linear form and consists of a list of Terms, each containing
    a pair of a ReferenceTensor and a GeometryTensor.

    Attributes:

        terms   - a list of Terms (products A0 * GS)
        aS      - a list of precomputed element tensor declarations
        a0      - a list of precomputed reference tensor declarations
        gS      - a list of precomputed geometry tensor declarations
        facet   - local number of facet
        num_ops - number of operations in computation of element tensor
    """

    def __init__(self, sum, format, cS_used, gS_used, options, facet):
        "Create ElementTensor."

        # Reset number of operations
        self.num_ops = 0

        # Compute terms
        self.terms = compute_terms(sum, Integral.EXTERIOR_FACET, facet, None, None)
        if len(self.terms) == 0:
            return

        # Compute element tensor declarations
        (self.aS, self.num_ops) = compute_element_tensor(self.terms, format, options)

        # Compute reference tensor declarations
        self.a0 = compute_reference_tensor(self.terms, format)
        
        # Compute geometry tensor declarations
        check_used(self.terms, format, gS_used)
        self.gS = compute_geometry_tensor(self.terms, format, gS_used, cS_used)

        # Save facet
        self.facet = facet
