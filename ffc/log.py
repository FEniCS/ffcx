"""This module provides functions used by the FFC implementation to
output messages. These may be redirected by the user of FFC.

This module reuses the corresponding log.py module from UFL which
is a wrapper for the standard Python logging module.
"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-01-12"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2009-12-17

# UFL modules.
from ufl.log import Logger
from ufl.log import log_functions
from ufl.log import DEBUG
from ufl.log import INFO
from ufl.log import ERROR

# Create FFC logger
ffc_logger = Logger("FFC")

# Create FFC global log functions
for foo in log_functions:
    exec("%s = lambda *message : ffc_logger.%s(*message)" % (foo, foo))

# Assertion, copied from UFL
def ffc_assert(condition, *message):
    "Assert that condition is true and otherwise issue an error with given message."
    condition or error(*message)

# Specialized FFC debugging tools

def debug_dict(d, title=""):
    "Pretty-print dictionary."
    if not title: title = "Dictionary"
    info("")
    begin(title)
    info("")
    for (key, value) in d.iteritems():
        info(key)
        info("-"*len(key))
        info(str(value))
        info("")
    end()

def debug_ir(ir, name=""):
    "Debug intermediate representation."
    title = "Intermediate representation"
    if name: title += " (%s)" % str(name)
    debug_dict(ir, title)

def debug_code(code, name=""):
    "Debug generated code."
    title = "Generated code"
    if name: title += " (%s)" % str(name)
    debug_dict(code, title)
