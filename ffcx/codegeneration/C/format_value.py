# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numbers
import re

_subs = (
    # Remove 0s after e+ or e-
    (re.compile(r"e[\+]0*(.)"), r"e\1"),
    (re.compile(r"e[\-]0*(.)"), r"e-\1"),
)


def format_float(x, precision=None):
    """Format a float value according to given precision."""
    global _subs

    if precision:
        if isinstance(x, complex):
            s = "({:.{prec}}+I*{:.{prec}})".format(x.real, x.imag, prec=precision)
        elif isinstance(x, float):
            s = "{:.{prec}}".format(x, prec=precision)
        else:
            s = "{:.{prec}}".format(float(x), prec=precision)
    else:
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
        raise RuntimeError("Unexpected type %s:\n%s" % (type(value),
                                                        str(value)))
