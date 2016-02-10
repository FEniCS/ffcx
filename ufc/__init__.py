__author__ = "Johan Hake (hake.dev@gmail.com)"
__copyright__ = "Copyright (C) 2010-2015 Johan Hake"
__date__ = "2010-08-19 -- 2015-02-26"
__license__  = "Released to the public domain"

# Modified by Anders Logg 2015

# Import Python versions of the abstract classes in the UFC interface
from .ufc import (cell,
                 function,
                 form,
                 finite_element,
                 dofmap,
                 cell_integral,
                 exterior_facet_integral,
                 interior_facet_integral,
                 vertex_integral,
                 custom_integral,
                 cutcell_integral,
                 interface_integral,
                 overlap_integral,
                 __version__,
                 __swigversion__)
