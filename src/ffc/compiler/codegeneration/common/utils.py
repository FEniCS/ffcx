"Code generation utils"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-04-11 -- 2007-04-11"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.constants import *

def inner_product(a, b, format):
    """Generate code for inner product of a and b, where a is a list
    of floating point numbers and b is a list of symbols."""

    # Check input
    if not len(a) == len(b):
        raise runtimeError, "Dimensions don't match for inner product."

    # Prefetch formats to speed up code generation
    format_add            = format["add"]
    format_subtract       = format["subtract"]
    format_multiply       = format["multiply"]
    format_floating_point = format["floating point"]
    
    # Add all entries
    value = None
    for i in range(len(a)):

        # Skip terms where a is almost zero
        if abs(a[i]) <= FFC_EPSILON:
            continue

        # Fancy handling of +, -, +1, -1
        if value:
            if abs(a[i] - 1.0) < FFC_EPSILON:
                value = format_add([value, b[i]])
            elif abs(a[i] + 1.0) < FFC_EPSILON:
                value = format_subtract([value, b[i]])
            elif a[i] > 0.0:
                value = format_add([value, format_multiply([format_floating_point(a[i]), b[i]])])
            else:
                value = format_subtract([value, format_multiply([format_floating_point(-a[i]), b[i]])])
        else:
            if abs(a[i] - 1.0) < FFC_EPSILON or abs(a[i] + 1.0) < FFC_EPSILON:
                value = b[i]
            else:
                value = format_multiply([format_floating_point(a[i]), b[i]])

    return value or format_floating_point(0.0)
