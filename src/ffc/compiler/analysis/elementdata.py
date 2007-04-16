__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-03-30 -- 2007-03-30"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *

# FFC fem modules
from ffc.fem.dofmap import *

class ElementData:
    """This class holds meta data for a list of elements. It has the
    same attributes as FormData, but only the elements and their dof
    maps are available."""

    def __init__(self, elements):
        "Create element for list of elements"

        debug("Extracting element data...")

        self.form                         = None
        self.signature                    = None
        self.rank                         = -1
        self.num_coefficients             = 0
        self.num_arguments                = len(elements)
        self.num_terms                    = 0
        self.num_cell_integrals           = 0
        self.num_exterior_facet_integrals = 0
        self.num_interior_facet_integrals = 0
        self.elements                     = elements
        self.dof_maps                     = [DofMap(element) for element in elements]
        self.cell_dimension               = 0

        debug("done")
