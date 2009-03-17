__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-09 -- 2009-03-09"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard and Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

def generate_reset_tensor(num_entries, format):
    "Generate code for resetting the entries of the local element tensor."

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

def generate_combined_code(codes, form_data, prefix, format):
    "Generate combined codes for tabulation of integrals."

    combined_code = {}

    # Generate combined code for cell integrals
    args = (codes, form_data, prefix, format, "cell_integral", form_data.num_cell_integrals)
    combined_code["cell_integrals"] = _generate_combined_code_common(*args)

    # Generate combined code for exterior facet integrals
    args = (codes, form_data, prefix, format, "exterior_facet_integral", form_data.num_exterior_facet_integrals)
    combined_code["exterior_facet_integrals"] = _generate_combined_code_common(*args)

    # Generate combined code for interior facet integrals
    args = (codes, form_data, prefix, format, "interior_facet_integral", form_data.num_interior_facet_integrals)
    combined_code["interior_facet_integrals"] = _generate_combined_code_common(*args)

    return combined_code

def _generate_combined_code_common(codes, form_data, prefix, format, integral_type, num_integrals):
    "Generate combined codes for tabulation of integrals."

    combined_code = []

    # Iterate over sub domains
    for i in range(num_integrals):

        # Add code for all representations
        contributions = []
        for (j, code) in enumerate(codes):
            key = (integral_type, i)
            if key in code:
                postfix = "%d_%s" % (i, code["representation"])
                combined_code.append((postfix, code[key]))
                contributions.append(postfix)

        # Add common code to sum up all contributions
        code = _generate_total_integral(integral_type, contributions, form_data, prefix, format)
        postfix = "%d" % i
        combined_code.append((postfix, code))

    return combined_code

def _generate_total_integral(integral_type, contributions, form_data, prefix, format):
    "Generate code for total tensor, summing contributions."

    code = {}

    # Add members
    code["members"] = ["\nprivate:\n"]
    for postfix in contributions:
        code["members"].append("  <form prefix>_%s_%s integral_%s;" % (integral_type, postfix, postfix))
    code["members"].append("")

    # FIXME: Possible optimization here not to reset entries if not needed

    # Reset all entries
    code["tabulate_tensor"] = []
    if integral_type == "interior_facet_integral":
        code["tabulate_tensor"] += generate_reset_tensor(form_data.num_entries_interior, format)
    else:
        code["tabulate_tensor"] += generate_reset_tensor(form_data.num_entries, format)

    # Sum contributions
    code["tabulate_tensor"].append("")
    code["tabulate_tensor"].append(format["comment"]("Add all contributions to element tensor"))
    for postfix in contributions:
        # FIXME: This is UFC specific
        if integral_type == "cell_integral":
            code["tabulate_tensor"] += ["integral_%s.tabulate_tensor(A, w, c);" % postfix]
        elif integral_type == "exterior_facet_integral":
            code["tabulate_tensor"] += ["integral_%s.tabulate_tensor(A, w, c, facet);" % postfix]
        else:
            code["tabulate_tensor"] += ["integral_%s.tabulate_tensor(A, w, c0, c1, facet0, facet1);" % postfix]

    return code
