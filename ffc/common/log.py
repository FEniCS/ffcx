"""This module provides functions used by the FFC implementation to
output messages. These may be redirected by the user of FFC.

This module reuses the corresponding log.py module from UFL which
is a wrapper for the standard Python logging module.
"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-01-12 -- 2009-03-06"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from ufl.log import *

# Create FFC logger
ffc_logger = Logger("FFC")

# Create FFC global log functions
for foo in log_functions:
    exec("%s = lambda *message : ffc_logger.%s(*message)" % (foo, foo))

# Assertion, copied from UFL
def ffc_assert(condition, *message):
    "Assert that condition is true and otherwise issue an error with given message."
    condition or error(*message)
