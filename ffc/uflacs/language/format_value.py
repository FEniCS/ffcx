# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
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

import numbers
import re

_subs = (
    # Remove 0s after e+ or e-
    (re.compile(r"e[\+]0*(.)"), r"e\1"),
    (re.compile(r"e[\-]0*(.)"), r"e-\1"),
    )
def format_float(x, precision=None):
    "Format a float value according to given precision."
    global _subs
    if precision:
        s = "{:.{prec}}".format(float(x), prec=precision)
    else:
        # Using "{}".format(float(x)) apparently results
        # in lesser precision in python 2 than python 3
        s = repr(float(x))
    for r, v in _subs:
        s = r.sub(v, s)
    return s


def format_int(x, precision=None):
    return str(x)


def format_value(value, precision=None):
    """Format a literal value as s tring.

    - float: Formatted according to current precision configuration.

    - int: Formatted as regular base 10 int literal.

    - str: Wrapped in "quotes".

    """
    if isinstance(value, numbers.Real):
        return format_float(float(value), precision=precision)
    elif isinstance(value, numbers.Integral):
        return format_int(int(value))
    elif isinstance(value, str):
        # FIXME: Is this ever used?
        assert '"' not in value
        return '"' + value + '"'
    elif hasattr(value, "ce_format"):
        return value.ce_format()
    else:
        raise RuntimeError("Unexpected type %s:\n%s" % (type(value), str(value)))
