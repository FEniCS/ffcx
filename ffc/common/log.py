"""This module provides functions used by the FFC implementation to
output messages. These may be redirected by the user of FFC.

This module reuses the corresponding log.py module from UFL which
is a wrapper for the standard Python logging module.
"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-01-12 -- 2009-01-12"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Import UFL log functions
from ufl.log import *

# Set up logger
set_logger("FFC")
