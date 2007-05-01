"Code generation for the UFC 1.0 format"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-08 -- 2007-03-27"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# UFC code templates
from ufc import *

# FFC common modules
from ffc.common.utils import *
from ffc.common.debug import *
from ffc.common.constants import *

# FFC language modules
from ffc.compiler.language.restriction import *
from ffc.compiler.language.integral import *

# FFC format modules
from codesnippets import *
from removeunused import *

# Choose map from restriction
choose_map = {Restriction.PLUS: "0", Restriction.MINUS: "1", None: ""}
# Transform format options  based on the sign of the power of the transform:
transform_options = {True: lambda m, j, k: "Jinv%s_%d%d" % (m, j, k),
                     False: lambda m, j, k: "J%s_%d%d" % (m, k, j)}
# Options for the printing q or 1.0/q for q string:
power_options = {True: lambda q: q, False: lambda q: "1.0/(%s)" % q}

# Specify formatting for code generation
format = { "add": lambda v: " + ".join(v),
           "subtract": lambda v: " - ".join(v),
           "multiply": lambda v: "*".join(v),
           "grouping": lambda v: "(%s)" % v,
           "block": lambda v: "{%s}" % v,
           "block begin": "{",
           "block end": "}",
           "separator": ", ",
           "return": lambda v: "return %s;" % v,
           "bool": lambda b: {True: "true", False: "false"}[b],
           "floating point": lambda a: "%.12g" % a,
           "tmp declaration": lambda j, k: "const double " + format["tmp access"](j, k),
           "tmp access": lambda j, k: "tmp%d_%d" % (j, k),
           "comment": lambda s: "// %s" % s,
           "determinant": "detJ",
           "scale factor": "det",
           "power": lambda base, exp: power_options[exp >= 0](format["multiply"]([str(base)]*abs(exp))),
           "constant": lambda j: "c%d" % j,
           "coefficient table": lambda j, k: "w[%d][%d]" % (j, k),
           "coefficient": lambda j, k: "w[%d][%d]" % (j, k),
           "transform": lambda p, j, k, r: "%s" % (transform_options[p >= 0](choose_map[r], j, k)),
           "reference tensor" : lambda j, i, a: None,
           "geometry tensor declaration": lambda j, a: "const double " + format["geometry tensor access"](j, a),
           "geometry tensor access": lambda j, a: "G%d_%s" % (j, "_".join(["%d" % index for index in a])),
           "element tensor": lambda i, k: "A[%d]" % k,
           "vertex values": lambda i: "vertex_values[%d]" % i,
           "dof values": lambda i: "dof_values[%d]" % i,
           "dofs": lambda i: "dofs[%d]" % i,
           "entity index": lambda d, i: "c.entity_indices[%d][%d]" % (d, i),
           "num entities": lambda dim : "m.num_entities[%d]" % dim,
           "offset declaration": "unsigned int offset",
           "offset access": "offset",
           "cell shape": lambda i: {1: "ufc::line", 2: "ufc::triangle", 3: "ufc::tetrahedron"}[i],
# quadrature, evalutate_basis(), evaluate_basis_derivatives()
# declarations
           "float declaration": "double ",
           "const float declaration": "const double ",
           "uint declaration": "unsigned int ",
           "const uint declaration": "const unsigned int ",
           "table declaration": "const static double ",
# access
           "array access": lambda i: "[%s]" %(i),
           "matrix access": lambda i,j: "[%s][%s]" %(i,j),
           "secondary index": lambda i: "_%d" %(i),
# program flow
           "dof map if": lambda i,j: "if (%d <= %s and %s <= %d)" %(i,\
                         format["argument basis num"], format["argument basis num"], j),
           "loop": lambda i,j,k: "for (unsigned int %s = %s; %s < %s; %s++)"% (i, j, i, k, i),
           "if": "if",
# operators
           "times equal": lambda i,j: "%s *= %s;" %(i,j),
           "add equal": lambda i,j: "%s += %s;" % (i,j),
           "is equal": " == ",
           "absolute value": lambda i: "std::abs(%s)" % (i),
# variable names
           "element tensor quad": lambda k: "A[%s]" % k,
           "loop integration points": "ip",
           "derivatives": lambda i,j,k,l: "dNdx%d_%d[%s][%s]" % (i,j,k,l),
           "element coordinates": lambda i,j: "x[%s][%s]" % (i,j),
           "weights": lambda i,j: "Weight%d[%s]" % (i,j),
           "psis": lambda l: "Psi%d_%d_%d_%s[%s][%s]" % (l[0],l[1],l[2],\
                             "".join(["%d" % index for index in l[3]]), l[4], l[5]),
           "argument coordinates": "coordinates",
           "argument values": "values",
           "argument basis num": "i",
           "argument derivative order": "n",
           "local dof": "dof",
           "x coordinate": "x",
           "y coordinate": "y",
           "z coordinate": "z",
           "scalings": lambda i,j: "scalings_%s_%d" %(i,j),
           "coefficients table": lambda i: "coefficients%d" %(i),
           "dmats table": lambda i: "dmats%d" %(i),
           "coefficient scalar": lambda i: "coeff%d" %(i),
           "new coefficient scalar": lambda i: "new_coeff%d" %(i),
           "psitilde_a": "psitilde_a",
           "psitilde_bs": lambda i: "psitilde_bs_%d" %(i),
           "psitilde_cs": lambda i,j: "psitilde_cs_%d%d" %(i,j),
           "basisvalue": lambda i: "basisvalue%d" %(i),
           "num derivatives": "num_derivatives",
           "reference derivatives": "derivatives",
           "derivative combinations": "combinations",
           "transform matrix": "transform",
           "transform Jinv": "Jinv",
# snippets
           "coordinate map": lambda i: {2:map_coordinates_2D, 3:map_coordinates_3D}[i],
           "snippet dof map": evaluate_basis_dof_map,
           "snippet eta_triangle": eta_triangle_snippet,
           "snippet eta_tetrahedron": eta_tetrahedron_snippet,
           "snippet combinations": combinations_snippet,
           "snippet transform2D": transform2D_snippet,
           "snippet transform3D": transform3D_snippet,
           "snippet evaluate_dof": lambda d : {2: evaluate_dof_2D, 3: evaluate_dof_3D}[d],
           "get cell vertices" : "const double * const * x = c.coordinates;",
# misc
           "block separator": ",\n",
           "new line": "\\\n",
           "end line": ";",
           "pointer": "*",
           "new": "new ",
           "delete": "delete "}

def write(generated_forms, prefix, options):
    "Generate UFC 1.0 code for a given list of pregenerated forms"
    debug("Generating code for UFC 1.0")

    # Generate code for header
    output = ""
    output += generate_header(prefix, options)
    output += "\n"

    # Generate UFC code
    output += generate_ufc(generated_forms, prefix, options)

    # Generate code for footer
    output += generate_footer(prefix, options)

    # Write file
    filename = "%s.h" % prefix
    file = open(filename, "w")
    file.write(output)
    file.close()
    debug("Output written to " + filename)

def generate_header(prefix, options):
    "Generate code for header"

    # Check if BLAS is required
    if options["blas"]:
        blas_include = "\n#include <cblas.h>"
        blas_warning = "\n// Warning: This code was generated with '-f blas' and requires cblas.h."
    else:
        blas_include = ""
        blas_warning = ""
        
    return """\
// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version %s.%s

#ifndef __%s_H
#define __%s_H

#include <cmath>
#include <ufc.h>%s
""" % (FFC_VERSION, blas_warning, prefix.upper(), prefix.upper(), blas_include)

def generate_footer(prefix, options):
    "Generate code for footer"
    return """\
#endif
"""

def generate_ufc(generated_forms, prefix, options):
    "Generate code for body"

    output = ""
    
    # Iterate over forms
    for i in range(len(generated_forms)):

        # Get pregenerated code, form data and prefix
        (form_code, form_data) = generated_forms[i]
        form_prefix = compute_prefix(prefix, generated_forms, i)

        # Generate code for ufc::finite_element(s)
        for (label, sub_element) in form_code["finite_elements"]:
            output += __generate_finite_element(sub_element, form_data, options, form_prefix, label)
            output += "\n"

        # Generate code for ufc::dof_map(s)
        for (label, sub_dof_map) in form_code["dof_maps"]:
            output += __generate_dof_map(sub_dof_map, form_data, options, form_prefix, label)
            output += "\n"

        # Generate code for ufc::cell_integral
        for j in range(form_data.num_cell_integrals):
            output += __generate_cell_integral(form_code[("cell_integral", j)], form_data, options, form_prefix, j)
            output += "\n"

        # Generate code for ufc::exterior_facet_integral
        for j in range(form_data.num_exterior_facet_integrals):
            output += __generate_exterior_facet_integral(form_code[("exterior_facet_integral", j)], form_data, options, form_prefix, j)
            output += "\n"
    
        # Generate code for ufc::interior_facet_integral
        for j in range(form_data.num_interior_facet_integrals):
            output += __generate_interior_facet_integral(form_code[("interior_facet_integral", j)], form_data, options, form_prefix, j)
            output += "\n"

        # Generate code for ufc::form
        if "form" in form_code:
            output += __generate_form(form_code["form"], form_data, options, form_prefix)
            output += "\n"

    return output

def compute_prefix(prefix, generated_forms, i):
    "Compute prefix for form i"

    # Get form ranks
    ranks = [form_data.rank for (form_code, form_data) in generated_forms]

    # Return prefixFunctional, prefixLinearForm or prefixBilinearForm
    # when we have exactly one form of ranks 0, 1 or 2
    count = [ranks.count(0), ranks.count(1), ranks.count(2)]
    if len(ranks) <= 3 and sum(count) > 0 and min(count) >= 0 and max(count) <= 1:
        postfixes = ["Functional", "LinearForm", "BilinearForm"]
        return "%s%s" % (prefix, postfixes[ranks[i]])

    # Return prefix_i if we have more than one rank
    if len(ranks) > 1:
        return "%s_%d" % (prefix, i)

    # Else, just return prefix
    return prefix

def __generate_finite_element(code, form_data, options, prefix, label):
    "Generate code for ufc::finite_element"

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = "%s_finite_element_%s" % (prefix, "_".join([str(i) for i in label]))

    # Generate code for members
    ufc_code["members"] = ""

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for signature
    ufc_code["signature"] = "return \"%s\";" % code["signature"]

    # Generate code for cell_shape
    ufc_code["cell_shape"] = "return %s;" % code["cell_shape"]
    
    # Generate code for space_dimension
    ufc_code["space_dimension"] = "return %s;" % code["space_dimension"]

    # Generate code for value_rank
    ufc_code["value_rank"] = "return %s;" % code["value_rank"]

    # Generate code for value_dimension
    cases = ["return %s;" % case for case in code["value_dimension"]]
    ufc_code["value_dimension"] = __generate_switch("i", cases, "return 0;")

    # Generate code for evaluate_basis
    ufc_code["evaluate_basis"] = __generate_body(code["evaluate_basis"])

    # Generate code for evaluate_basis
    ufc_code["evaluate_basis_derivatives"] = __generate_body(code["evaluate_basis_derivatives"])

    # Generate code for evaluate_dof
    ufc_code["evaluate_dof"] = __generate_body(code["evaluate_dof"])

    # Generate code for inperpolate_vertex_values
    ufc_code["interpolate_vertex_values"] = __generate_body(code["interpolate_vertex_values"])

    # Generate code for num_sub_elements
    ufc_code["num_sub_elements"] = "return %s;" % code["num_sub_elements"]

    # Generate code for sub_element
    num_sub_elements = eval(code["num_sub_elements"])
    if num_sub_elements == 1:
        ufc_code["create_sub_element"] = "return new %s();" % ufc_code["classname"]
    else:
        cases = ["return new %s_%d();" % (ufc_code["classname"], i) for i in range(num_sub_elements)]
        ufc_code["create_sub_element"] = __generate_switch("i", cases, "return 0;")

    return __generate_code(finite_element_combined, ufc_code)

def __generate_dof_map(code, form_data, options, prefix, label):
    "Generate code for ufc::dof_map"

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = "%s_dof_map_%s" % (prefix, "_".join([str(i) for i in label]))

    # Generate code for members
    ufc_code["members"] = "\nprivate:\n\n  unsigned int __global_dimension;\n"

    # Generate code for constructor
    ufc_code["constructor"] = "__global_dimension = 0;"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for signature
    ufc_code["signature"] = "return \"%s\";" % code["signature"]

    # Generate code for needs_mesh_entities
    cases = ["return %s;" % case for case in code["needs_mesh_entities"]]
    ufc_code["needs_mesh_entities"] = __generate_switch("d", cases, "return false;")

    # Generate code for init_mesh
    ufc_code["init_mesh"] = "__global_dimension = %s;\nreturn false;" % code["global_dimension"]

    # Generate code for init_cell
    ufc_code["init_cell"] = "// Do nothing"

    # Generate code for init_cell_finalize
    ufc_code["init_cell_finalize"] = "// Do nothing"

    # Generate code for global_dimension
    ufc_code["global_dimension"] = "return __global_dimension;"

    # Generate code for local dimension
    ufc_code["local_dimension"] = "return %s;" % code["local_dimension"]

    # Generate code for num_facet_dofs
    ufc_code["num_facet_dofs"] = "return %s;" % code["num_facet_dofs"]

    # Generate code for tabulate_dofs
    ufc_code["tabulate_dofs"] = __generate_body(code["tabulate_dofs"])

    # Generate code for tabulate_facet_dofs
    ufc_code["tabulate_facet_dofs"] = __generate_switch("facet", [__generate_body(case) for case in code["tabulate_facet_dofs"]])

    # Generate code for tabulate_coordinates
    ufc_code["tabulate_coordinates"] = "// Not implemented"
    #ufc_code["tabulate_coordinates"] = __generate_body(code["tabulate_coordinates"])

    # Generate code for num_sub_dof_maps
    ufc_code["num_sub_dof_maps"] = "return %s;" % code["num_sub_dof_maps"]

    # Generate code for create_sub_dof_map
    num_sub_dof_maps = eval(code["num_sub_dof_maps"])
    if num_sub_dof_maps == 1:
        ufc_code["create_sub_dof_map"] = "return new %s();" % ufc_code["classname"]
    else:
        cases = ["return new %s_%d();" % (ufc_code["classname"], i) for i in range(num_sub_dof_maps)]
        ufc_code["create_sub_dof_map"] = __generate_switch("i", cases, "return 0;")


    return __generate_code(dof_map_combined, ufc_code)

def __generate_cell_integral(code, form_data, options, prefix, i):
    "Generate code for ufc::cell_integral"

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = "%s_cell_integral_%d" % (prefix, i)

    # Generate code for members
    ufc_code["members"] = ""

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for tabulate_tensor
    body  = __generate_jacobian(form_data.cell_dimension, Integral.CELL)
    body += "\n"
    body += __generate_body(code["tabulate_tensor"])
    ufc_code["tabulate_tensor"] = remove_unused(body)

    return __generate_code(cell_integral_combined, ufc_code)

def __generate_exterior_facet_integral(code, form_data, options, prefix, i):
    "Generate code for ufc::exterior_facet_integral"

    ufc_code = {}
    
    # Set class name
    ufc_code["classname"] = "%s_exterior_facet_integral_%d" % (prefix, i)

    # Generate code for members
    ufc_code["members"] = ""

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for tabulate_tensor
    switch = __generate_switch("facet", [__generate_body(case) for case in code["tabulate_tensor"][1]])
    body  = __generate_jacobian(form_data.cell_dimension, Integral.EXTERIOR_FACET)
    body += "\n"
    body += __generate_body(code["tabulate_tensor"][0])
    body += "\n"
    body += switch
    ufc_code["tabulate_tensor"] = remove_unused(body)
    
    return __generate_code(exterior_facet_integral_combined, ufc_code)

def __generate_interior_facet_integral(code, form_data, options, prefix, i):
    "Generate code for ufc::interior_facet_integral"

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = "%s_interior_facet_integral_%d" % (prefix, i)

    # Generate code for members
    ufc_code["members"] = ""

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for tabulate_tensor, impressive line of Python code follows
    switch = __generate_switch("facet0", [__generate_switch("facet1", [__generate_body(case) for case in cases]) for cases in code["tabulate_tensor"][1]])
    body  = __generate_jacobian(form_data.cell_dimension, Integral.INTERIOR_FACET)
    body += "\n"
    body += __generate_body(code["tabulate_tensor"][0])
    body += "\n"
    body += switch
    ufc_code["tabulate_tensor"] = remove_unused(body)

    return __generate_code(interior_facet_integral_combined, ufc_code)

def __generate_form(code, form_data, options, prefix):
    "Generate code for ufc::form"

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = prefix

    # Generate code for members
    ufc_code["members"] = ""

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for signature
    ufc_code["signature"] = "return \"%s\";" % code["signature"]

    # Generate code for rank
    ufc_code["rank"] = "return %s;" % code["rank"]

    # Generate code for num_coefficients
    ufc_code["num_coefficients"] = "return %s;" % code["num_coefficients"]
    
    # Generate code for num_cell_integrals
    ufc_code["num_cell_integrals"] = "return %s;" % code["num_cell_integrals"]

    # Generate code for num_exterior_facet_integrals
    ufc_code["num_exterior_facet_integrals"] = "return %s;" % code["num_exterior_facet_integrals"]
    
    # Generate code for num_interior_facet_integrals
    ufc_code["num_interior_facet_integrals"] = "return %s;" % code["num_interior_facet_integrals"]

    # Generate code for create_finite_element
    num_cases = form_data.num_arguments
    cases = ["return new %s_finite_element_%d();" % (prefix, i) for i in range(num_cases)]
    ufc_code["create_finite_element"] = __generate_switch("i", cases, "return 0;")

    # Generate code for create_dof_map
    num_cases = form_data.num_arguments
    cases = ["return new %s_dof_map_%d();" % (prefix, i) for i in range(num_cases)]
    ufc_code["create_dof_map"] = __generate_switch("i", cases, "return 0;")

    # Generate code for cell_integral
    num_cases = form_data.num_cell_integrals
    cases = ["return new %s_cell_integral_%d();" % (prefix, i) for i in range(num_cases)]
    ufc_code["create_cell_integral"] = __generate_switch("i", cases, "return 0;")

    # Generate code for exterior_facet_integral
    num_cases = form_data.num_exterior_facet_integrals
    cases = ["return new %s_exterior_facet_integral_%d();" % (prefix, i) for i in range(num_cases)]
    ufc_code["create_exterior_facet_integral"] = __generate_switch("i", cases, "return 0;")

    # Generate code for interior_facet_integral
    num_cases = form_data.num_interior_facet_integrals
    cases = ["return new %s_interior_facet_integral_%d();" % (prefix, i) for i in range(num_cases)]
    ufc_code["create_interior_facet_integral"] = __generate_switch("i", cases, "return 0;")

    return __generate_code(form_combined, ufc_code)

def __generate_jacobian(cell_dimension, integral_type):
    "Generate code for computing jacobian"

    # Choose space dimension
    if cell_dimension == 2:
        jacobian = jacobian_2D
        facet_determinant = facet_determinant_2D
    else:
        jacobian = jacobian_3D
        facet_determinant = facet_determinant_2D
    
    # Check if we need to compute more than one Jacobian
    if integral_type == Integral.CELL:
        code  = jacobian % {"restriction":  ""}
        code += "\n"
        code += scale_factor
    elif integral_type == Integral.EXTERIOR_FACET:
        code  = jacobian % {"restriction":  ""}
        code += "\n"
        code += facet_determinant % {"restriction": "", "facet" : "facet"}
    elif integral_type == Integral.INTERIOR_FACET:
        code  = jacobian % {"restriction": choose_map[Restriction.PLUS]}
        code += "\n"
        code += jacobian % {"restriction": choose_map[Restriction.MINUS]}
        code += "\n"
        code += facet_determinant % {"restriction": choose_map[Restriction.PLUS], "facet": "facet0"}

    return code

def __generate_switch(variable, cases, default = ""):
    "Generate switch statement from given variable and cases"

    # Special case: no cases
    if len(cases) == 0:
        return default

    # Special case: one case
    if len(cases) == 1:
        return cases[0]

    # Create switch
    code = "switch ( %s )\n{\n" % variable
    for i in range(len(cases)):
        code += "case %d:\n%s\n  break;\n" % (i, indent(cases[i], 2))
    code += "}"
    if not default == "":
        code += "\n" + default
    
    return code

def __generate_body(declarations):
    "Generate function body from list of declarations or statements"
    lines = []
    for declaration in declarations:
        if isinstance(declaration, tuple):
            lines += ["%s = %s;" % declaration]
        else:
            lines += ["%s" % declaration]
    return "\n".join(lines)

def __generate_code(format_string, code):
    "Generate code according to format string and code dictionary"

    # Fix indentation
    for key in code:
        if not key in ["classname", "members"]:
            code[key] = indent(code[key], 4)

    # Generate code
    return format_string % code
