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

    # Generate code for cell integrals
    combined_code["cell_integrals"] = []
    for i in range(form_data.num_cell_integrals):

        # Add code for all representations
        contributions = []
        for (j, code) in enumerate(codes):
            key = ("cell_integral", i)
            if key in code:
                postfix = "%d_%s" % (i, code["representation"])
                combined_code["cell_integrals"].append((postfix, code[key]))
                contributions.append(postfix)

        # Add common code to sum up all contributions
        code = _generate_total_integral("cell_integral",
                                        contributions,
                                        form_data.num_entries,
                                        prefix,
                                        format)
        combined_code["cell_integrals"].append(("%d" % i, code))

    # Generate code for exterior facet integrals
    combined_code["exterior_facet_integrals_integrals"] = []
    for i in range(form_data.num_exterior_facet_integrals):

        # Add code for all representations
        for (j, code) in enumerate(codes):
            key = ("exterior_facet_integral", i)
            if key in code:
                postfix = "%d_%s" % (i, code["representation"])
                combined_code["exterior_facet_integrals"].append((postfix, code[key]))
                contributions.append(postfix)

        # Add common code to sum up all contributions
        code = _generate_total_integral("exterior_facet_integral",
                                        contributions,
                                        form_data.num_entries,
                                        prefix,
                                        format)
        combined_code["cell_integrals"].append(("%d" % i, code))

    # Generate code for interior facet integrals
    combined_code["interior_facet_integrals_integrals"] = []
    for i in range(form_data.num_interior_facet_integrals):

        # Add code for all representations
        for (j, code) in enumerate(codes):
            key = ("interior_facet_integral", i)
            if key in code:
                postfix = "%d_%s" % (i, code["representation"])
                combined_code["interior_facet_integrals"].append((postfix, code[key]))
                contributions.append(postfix)

        # Add common code to sum up all contributions
        code = _generate_total_integral("interior_facet_integral",
                                        contributions,
                                        form_data.num_entries,
                                        prefix,
                                        format)
        combined_code["interior_facet_integral"].append(("%d" % i, code))
       
    return combined_code

def _generate_total_integral(integral_type, contributions, num_entries, prefix, format):
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
    code["tabulate_tensor"] += generate_reset_tensor(num_entries, format)

    # Sum contributions
    code["tabulate_tensor"].append("")
    code["tabulate_tensor"].append(format["comment"]("Add all contributions to element tensor"))
    for postfix in contributions:
        code["tabulate_tensor"] += ["integral_%s.tabulate_tensor(A, w, c);" % postfix]

    return code
