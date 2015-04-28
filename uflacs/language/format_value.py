# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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

import re

_float_threshold = None
_float_precision = None
_float_fmt = None

def set_float_precision(precision, threshold=None):
    "Configure float formatting precision and zero threshold."
    global _float_epsilon, _float_precision, _float_fmt
    _float_precision = precision
    _float_threshold = threshold
    _float_fmt = "{{:.{0:d}e}}".format(_float_precision)

def reset_float_precision():
    "Set float precision and zero threshold back to default."
    set_float_precision(15, 1e-15)

# Execute default on startup
reset_float_precision()

_p0 = re.compile("0+e")
_p1 = re.compile("e\\+00$")
_p2 = re.compile("\\.$")
def format_float(x):
    "Format a float value according to set_float_precision."
    global _float_threshold, _float_fmt, _p0, _p1, _p2
    if _float_threshold is not None and abs(x) < _float_threshold:
        return "0.0"
    else:
        s = _float_fmt.format(x)
        s = s.strip()
        s = _p0.sub("e", s)
        s = _p1.sub("", s)
        s = _p2.sub(".0", s)
        return s

def format_value(snippets):
    """Format a simple value.

    - int or float: Formatted according to current configuration.

    - str: Used directly.

    """
    if isinstance(snippets, float):
        return format_float(snippets)
    elif isinstance(snippets, int):
        return str(snippets)
    elif isinstance(snippets, str):
        return snippets
    else:
        raise RuntimeError("Unexpected type %s:\n%s" % (type(snippets), str(snippets)))
