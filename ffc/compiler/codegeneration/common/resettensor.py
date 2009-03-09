__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-09 -- 2009-03-09"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard and Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

def generate_reset_tensor(num_entries, format):
    "Generate code for resetting the entries of the local element tensor."

    print "HEJ"

    # Generate code as a list of declarations
    code = []    

    # Comment
    code.append(format["comment"]("Reset values of the element tensor block"))

    # Prefetch formats to speed up code generation
    format_element_tensor = format["element tensor"]
    format_floating_point = format["floating point"]

    # Set entries to zero
    for k in range(num_entries):
        name = format_element_tensor(None, k)
        value = format_floating_point(0.0)
        code += [(name, value)]

    return code
