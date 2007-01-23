__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2006-12-01"
__copyright__ = "Copyright (C) 2006 Anders Logg"
__license__  = "GNU GPL Version 2"

import integral

# Predefined integrals
dx = integral.Integral("cell")
ds = integral.Integral("exterior facet")
dS = integral.Integral("interior facet")
