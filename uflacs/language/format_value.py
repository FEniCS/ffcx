# -*- coding: utf-8 -*-
# Copyright (C) 2011-2016 Martin Sandve Alnæs
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
import numpy
#from ffc.log import info

_float_threshold = None
_float_precision = None
_float_fmt = "%r"


def set_float_precision(precision, threshold=None):
    "Configure float formatting precision and zero threshold."
    global _float_precision, _float_threshold, _float_fmt
    if threshold is None:
        threshold = 10.0**-(precision-1)  # Matching FFC behaviour, I'd like to drop the -1 here
    _float_precision = precision
    _float_threshold = threshold
    #_float_fmt = "{{:.{0:d}e}}".format(_float_precision)
    if _float_precision is None:
        _float_fmt = "%r"
    else:
        _float_fmt = "%%.%dg" % _float_precision
    #info("Setting float precision to %d in uflacs." % (precision,))


def get_float_precision():
    return _float_precision


def get_float_threshold():
    return _float_threshold


def reset_float_precision():
    "Set float precision and zero threshold back to default."
    set_float_precision(15)


# Execute default on startup
reset_float_precision()


_subs = (
    # Remove trailing spaces (using .strip() for this)
    #(re.compile(r"^[ ]*([^ ])[ ]*$"), "\1"),
    # Remove 0s before e
    (re.compile(r"0+e"), "e"),
    # Remove 0s after e+ or e-
    (re.compile(r"e([\+\-])0+"), r"e\1"),
    # Remove e+/e- at end (zero removal above can produce this)
    (re.compile(r"e[\+\-]$"), ""),
    # Replace lonely . at end with .0 for consistency
    (re.compile("\\.$"), ".0"),
    # Append .0 if it's an integer
    (re.compile(r"^([\+\-]{,1}\d+)$"), r"\1.0"),
    )
def format_float(x):
    "Format a float value according to set_float_precision."
    global _float_threshold, _float_fmt, _subs

    if _float_threshold is not None and abs(x) < _float_threshold:
        return "0.0"
    else:
        s = (_float_fmt % x).strip()
        for r, v in _subs:
            s = r.sub(v, s)
        return s


_ints = (int, numpy.integer)
_floats = (float, numpy.floating)
def format_value(value):
    """Format a literal value as s tring.

    - float: Formatted according to current precision configuration.

    - int: Formatted as regular base 10 int literal.

    - str: Wrapped in "quotes".

    """
    global _floats, _ints
    if isinstance(value, _floats):
        return format_float(float(value))
    elif isinstance(value, _ints):
        return str(int(value))
    elif isinstance(value, str):
        return '"' + value + '"'
    elif hasattr(value, "ce_format"):
        return value.ce_format()
    else:
        raise RuntimeError("Unexpected type %s:\n%s" % (type(value), str(value)))
