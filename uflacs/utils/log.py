"""This module provides functions used by the UFLACS implementation to
output messages. These may be redirected by the user of UFLACS.

This module reuses the corresponding log.py module from UFL which
is a wrapper for the standard Python logging module.
"""

# Copyright (C) 2009 Anders Logg
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Kristian B. Oelgaard, 2009
# Modified by Martin S. Alnaes, 2013
#
# First added:  2009-01-12
# Last changed: 2013-08-18

# UFL modules
from ufl.log import Logger
from ufl.log import log_functions
from ufl.log import INFO, DEBUG, ERROR, CRITICAL
from ufl.common import dstr, tstr

# Base class for UFLACS exceptions
class UFLACSException(Exception):
    "Base class for UFLACS exceptions"
    pass

# Create UFLACS logger
uflacs_logger = Logger("UFLACS", UFLACSException)

# Create UFLACS global log functions
for foo in log_functions:
    exec("%s = lambda *message : uflacs_logger.%s(*message)" % (foo, foo))

# Assertion, copied from UFL
def uflacs_assert(condition, *message):
    "Assert that condition is true and otherwise issue an error with given message."
    condition or error(*message)

# Set default log level
set_level(INFO)

