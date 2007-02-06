__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2007-02-06"
__copyright__ = "Copyright (C) 2006-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC language modules
from integral import *
from index import *

# Predefined integrals
dx = Integral("cell")
ds = Integral("exterior facet")
dS = Integral("interior facet")

# Predefined indices
i = Index()
j = Index()
k = Index()
l = Index()
m = Index()
n = Index()
