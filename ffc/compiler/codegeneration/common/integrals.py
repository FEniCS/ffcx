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

def combine_tensors(code, code0, code1, key, reset_code, facet_integral):
    "Combine code for tabulate_tensor from two different code generators."

    # If subdomain has both representations then combine them
    if key in code1 and key in code0:
        code[key] = {("tabulate_tensor_tensor"):code1[key],
                     ("tabulate_tensor_quadrature"):code0[key],
                     "reset_tensor": reset_code}

    # Handle code from tensor generator
    elif key in code1:
        # Add reset code to tabulate_tensor code
        val = code1[key]["tabulate_tensor"]
        # Check if we have a tuple (common, cases) for facet integrals
        if facet_integral:
            val = (reset_code + val[0], val[1])
        else:
            val = reset_code + val
        code1[key]["tabulate_tensor"] = val
        code[key] = code1[key]

    # Handle code from quadrature generator
    elif key in code0:
        # Add reset code to tabulate_tensor code
        val = code0[key]["tabulate_tensor"]
        # Check if we have a tuple (common, cases) for facet integrals
        if facet_integral:
            val = (reset_code + val[0], val[1])
        else:
            val = reset_code + val
        code0[key]["tabulate_tensor"] = val
        code[key] = code0[key]

    # If we reach this level it means that no code has been generated
    # for the given subdomain so we need to add the reset code
    # NOTE: If we were sure that all assemblers would reset the local
    # tensor before calling tabulate_tensor this wouldn't be needed
    else:
        if facet_integral:
            code[key] = {"tabulate_tensor": (reset_code, []), "members": []}
        else:
            code[key] = {"tabulate_tensor": reset_code, "members": []}
