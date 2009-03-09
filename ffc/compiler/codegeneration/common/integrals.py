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

def combine_tensors(code, quadrature_code, tensor_code, key, reset_code, facet_integral):
    "Combine code for tabulate_tensor for the two different code generators."

    # Note: This could be simplified if the assembler would reset
    # the tensor before calling tabulate_tensor

    # If subdomain has both representations then combine them
    if key in tensor_code and key in quadrature_code:
        code[key] = {("tabulate_tensor_quadrature"): quadrature_code[key],
                     ("tabulate_tensor_tensor"): tensor_code[key],
                     "reset_tensor": reset_code}

    # Only quadrature code generated
    elif key in quadrature_code:
        
        # Add reset code to tabulate_tensor code
        value = quadrature_code[key]["tabulate_tensor"]
        
        # Handle facet integral cases
        if facet_integral:
            value = (reset_code + value[0], value[1])
        else:
            value = reset_code + value

        quadrature_code[key]["tabulate_tensor"] = value
        code[key] = quadrature_code[key]

    # Only tensor code generated
    elif key in tensor_code:
        
        # Add reset code to tabulate_tensor code
        value = tensor_code[key]["tabulate_tensor"]

        # Handle facet integral cases
        if facet_integral:
            value = (reset_code + value[0], value[1])
        else:
            value = reset_code + value
            
        tensor_code[key]["tabulate_tensor"] = value
        code[key] = tensor_code[key]

    # No code generated
    else:
        if facet_integral:
            code[key] = {"tabulate_tensor": (reset_code, []), "members": []}
        else:
            code[key] = {"tabulate_tensor": reset_code, "members": []}
