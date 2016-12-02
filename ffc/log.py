# -*- coding: utf-8 -*-
"""This module provides functions used by the FFC implementation to
output messages. These may be redirected by the user of FFC.

This module reuses the corresponding log.py module from UFL which
is a wrapper for the standard Python logging module.
"""

# Copyright (C) 2009 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Kristian B. Oelgaard, 2009

# UFL modules
from ufl.log import Logger
from ufl.log import log_functions
from ufl.log import INFO, DEBUG, WARNING, ERROR, CRITICAL
from ufl.utils.sorting import sorted_by_key
from ufl.utils.formatting import dstr, tstr

# Create FFC logger
ffc_logger = Logger("FFC")

# Create FFC global log functions
for foo in log_functions:
    exec("%s = lambda *message : ffc_logger.%s(*message)" % (foo, foo))


# Assertion, copied from UFL
def ffc_assert(condition, *message):
    "Assert that condition is true and otherwise issue an error with given message."
    condition or error(*message)


# Set default log level
set_level(INFO)

#--- Specialized FFC debugging tools ---


def debug_dict(d, title=""):
    "Pretty-print dictionary."
    if not title:
        title = "Dictionary"
    info("")
    begin(title)
    info("")
    for (key, value) in sorted_by_key(d):
        info(key)
        info("-" * len(key))
        info(str(value))
        info("")
    end()


def debug_ir(ir, name=""):
    "Debug intermediate representation."
    title = "Intermediate representation"
    if name:
        title += " (%s)" % str(name)
    debug_dict(ir, title)


def debug_code(code, name=""):
    "Debug generated code."
    title = "Generated code"
    if name:
        title += " (%s)" % str(name)
    debug_dict(code, title)
