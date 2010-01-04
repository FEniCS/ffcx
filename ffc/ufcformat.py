"Code generation for the UFC 1.0 format"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-08"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.
# Modified by Dag Lindbo, 2008.
# Modified by Johan Hake, 2009.
# Modified by Garth N. Wells, 2009.
# Last changed: 2010-01-04

# Python modules.
import os
import platform

# UFC code templates
# UFC build
from ufc_utils import build_ufc_module

# UFC function.
from ufc_utils import function_combined
from ufc_utils import function_header
from ufc_utils import function_implementation

# UFC finite_element.
from ufc_utils import finite_element_combined
from ufc_utils import finite_element_header
from ufc_utils import finite_element_implementation

# UFC dof_map.
from ufc_utils import dof_map_combined
from ufc_utils import dof_map_header
from ufc_utils import dof_map_implementation

# UFC integrals.
from ufc_utils import cell_integral_combined
from ufc_utils import cell_integral_header
from ufc_utils import cell_integral_implementation
from ufc_utils import exterior_facet_integral_combined
from ufc_utils import exterior_facet_integral_header
from ufc_utils import exterior_facet_integral_implementation
from ufc_utils import interior_facet_integral_combined
from ufc_utils import interior_facet_integral_header
from ufc_utils import interior_facet_integral_implementation

# UFC form.
from ufc_utils import form_combined
from ufc_utils import form_header
from ufc_utils import form_implementation

# DOLFIN wrapper generator.
try:
    from dolfin_utils.wrappers import generate_dolfin_code
    from dolfin_utils.wrappers import UFCFormNames
    dolfin_utils_imported = True
except:
    dolfin_utils_imported = False

# FFC modules.
from log import info
from log import error
from constants import FFC_VERSION
from codegeneratorsutils import indent
from removeunused import remove_unused
from compiler import ElementData

# FFC format modules
# TODO: Finish this import list
#from codesnippets import evaluate_basis_dof_map, eta_interval_snippet, eta_triangle_snippet, eta_tetrahedron_snippet, combinations_snippet, calculate_dof, map_coordinates_interval, map_coordinates_triangle, map_coordinates_tetrahedron, transform_interval_snippet, transform_triangle_snippet, transform_tetrahedron_snippet, map_onto_physical_1D, map_onto_physical_2D, map_onto_physical_3D, jacobian_1D, jacobian_2D, jacobian_3D, facet_determinant_1D, facet_determinant_2D, facet_determinant_3D, scale_factor, normal_direction_2D

from codesnippets import *

# Choose map from restriction
choose_map = {None: "", "+": "0", "-": 1}
transform_options = {"JINV": lambda m, j, k: "Jinv%s_%d%d" % (m, j, k),
                     "J": lambda m, j, k: "J%s_%d%d" % (m, k, j)}
transform_options_ufl = {"JINV": lambda m, j, k: "Jinv%s_%d%d" % (m, j, k),
                         "J": lambda m, j, k: "J%s_%d%d" % (m, j, k)}
# Options for the printing q or 1.0/(q) for q string:
power_options = {True: lambda q: q, False: lambda q: "1.0/(%s)" % q}

class Format:

    def __init__(self, options):
        "Initialize code generation for given options"

        # Check format option
        self.output_format = options["format"].lower()
        if not self.output_format in ["ufc", "dolfin"]:
            error("Don't know how to compile code for format '%s'." % output_format)

        # Check that DOLFIN wrapper utils have been imorted
        if self.output_format == "dolfin" and not dolfin_utils_imported:
            error("Module dolfin_utils must be imported to generate wrapper for DOLFIN.")

        # Attach format
        # FIXME: KBO: It should be possible to clean up the below formats by only
        # having e.g., format["const"] = "const ", format["double"] = "double "
        # instead of format["const float declaration"] = "const double ".
        # Then when we generating code, we have to combine the two manually, it
        # requires a little more work when using the formt, but the format itself
        # will be a lot cleaner and less confusion is likely to appear.

    def write(self, generated_forms, prefix, options):
        "Generate UFC 1.0 code for a given list of pregenerated forms"

        # Strip directory names from prefix and add output directory
        prefix = prefix.split(os.path.join(' ',' ').split()[0])[-1]
        full_prefix = os.path.join(options["output_dir"], prefix)

        # Generate code for header
        output = ""
        if self.output_format == "ufc":
            output += _generate_header(prefix, options)
        elif self.output_format == "dolfin":
            output += _generate_dolfin_header(prefix, options)
        output += "\n"

        if not options["split"]:

            if self.output_format == "dolfin":

                # Generate UFC code
                output += _generate_ufc(generated_forms, prefix, options, "combined", self.format)

                # Generate code for DOLFIN wrappers
                output += _generate_dolfin_wrappers(generated_forms, prefix, options, self.format)

            elif self.output_format == "ufc":
                # Generate UFC code
                output += _generate_ufc(generated_forms, prefix, options, "combined", self.format)

            # Generate code for footer
            output += _generate_footer(prefix, options)

            # Write file
            filename = "%s.h" % full_prefix
            file = open(filename, "w")
            file.write(output)
            file.close()
            info("Output written to " + filename + ".")

        else:

            if self.output_format == "dolfin":

                # Generate UFC header code
                output += _generate_ufc(generated_forms, prefix, options, "header", self.format)

                # Generate code for DOLFIN wrappers
                output += _generate_dolfin_wrappers(generated_forms, prefix, options, self.format)

            elif self.output_format == "ufc":
                # Generate UFC code
                output += _generate_ufc(generated_forms, prefix, options, "header", self.format)

            # Generate code for footer
            output += _generate_footer(prefix, options)

            # Write file
            filename = "%s.h" % full_prefix
            file = open(filename, "w")
            file.write(output)
            file.close()
            info("Output written to " + filename + ".")

            output = ""

            # Generate UFC implementation code
            output += "#include \"%s.h\"\n" % prefix

            if self.output_format == "dolfin":
                output += _generate_ufc(generated_forms, prefix, options, "implementation", self.format)

            elif self.output_format == "ufc":
                output += _generate_ufc(generated_forms, prefix, options, "implementation", self.format)

            # Write file
            filename = "%s.cpp" % full_prefix
            file = open(filename, "w")
            file.write(output)
            file.close()
            info("Output written to " + filename + ".")

def _generate_header(prefix, options):
    "Generate code for header"

    return """\
// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version %s.

#ifndef __%s_H
#define __%s_H

#include <cmath>
#include <stdexcept>
#include <ufc.h>
""" % (FFC_VERSION, prefix.upper(), prefix.upper())

def _generate_dolfin_header(prefix, options):
    "Generate DOLFIN file header"

    return """\
// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version %s.
//
// Warning: This code was generated with the option '-l dolfin'
// and contains DOLFIN-specific wrappers that depend on DOLFIN.

#ifndef __%s_H
#define __%s_H

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <ufc.h>
    """ % (FFC_VERSION, prefix.upper(), prefix.upper())

def _generate_footer(prefix, options):
    "Generate code for footer"
    return """\
#endif
"""

def _generate_ufc(generated_forms, prefix, options, code_section, format):
    "Generate code for body"

    info("Formatting code for UFC 1.0.")

    output = ""

    # Iterate over forms
    for (i, (form_code, form_data)) in enumerate(generated_forms):

        # Generate code for ufc::finite_element(s)
        for (label, element) in form_code["finite_elements"]:
            output += _generate_finite_element(element, form_data, options, prefix, i, label, code_section, format)
            output += "\n"

        # Generate code for ufc::dof_map(s)
        for (label, dof_map) in form_code["dof_maps"]:
            output += _generate_dof_map(dof_map, form_data, options, prefix, i, label, code_section, format)
            output += "\n"

        # Generate code for ufc::cell_integral
        if form_code.has_key("cell_integrals"):
            for (label, code) in form_code["cell_integrals"]:
                output += _generate_cell_integral(code, form_data, options, prefix, i, label, code_section, format)
                output += "\n"

        # Generate code for ufc::exterior_facet_integral
        if form_code.has_key("exterior_facet_integrals"):
            for (label, code) in form_code["exterior_facet_integrals"]:
                output += _generate_exterior_facet_integral(code, form_data, options, prefix, i, label, code_section, format)
                output += "\n"

        # Generate code for ufc::interior_facet_integral
        if form_code.has_key("interior_facet_integrals"):
            for (label, code) in form_code["interior_facet_integrals"]:
                output += _generate_interior_facet_integral(code, form_data, options, prefix, i, label, code_section, format)
                output += "\n"

        # Generate code for ufc::form
        if "form" in form_code:
            output += _generate_form(form_code["form"], form_data, options, prefix, i, code_section, format)
            output += "\n"

    return output

def _generate_finite_element(code, form_data, options, prefix, i, label, code_section, format):
    "Generate code for ufc::finite_element"

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = format["classname finite_element"](prefix, i, label)

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
    ufc_code["value_dimension"] = _generate_switch("i", cases, "return 0;")

    # Generate code for evaluate_basis (and vectorised counterpart)
    ufc_code["evaluate_basis"] = _generate_body(code["evaluate_basis"])
    ufc_code["evaluate_basis_all"] = _generate_body(code["evaluate_basis_all"])

    # Generate code for evaluate_basis_derivatives (and vectorised counterpart)
    ufc_code["evaluate_basis_derivatives"] = _generate_body(code["evaluate_basis_derivatives"])
    ufc_code["evaluate_basis_derivatives_all"] = _generate_body(code["evaluate_basis_derivatives_all"])

    # Generate code for evaluate_dof
    ufc_code["evaluate_dof"] = _generate_body(code["evaluate_dof"])

    # Generate code for evaluate_dofs (introduced in UFC 1.1)
    ufc_code["evaluate_dofs"] = format["exception"]("Not implemented (introduced in UFC v1.1).")

    # Generate code for inperpolate_vertex_values
    ufc_code["interpolate_vertex_values"] = remove_unused(_generate_body(code["interpolate_vertex_values"]))

    # Generate code for num_sub_elements
    ufc_code["num_sub_elements"] = "return %s;" % code["num_sub_elements"]

    # Generate code for sub_element
    num_sub_elements = eval(code["num_sub_elements"])
    if num_sub_elements == 1:
        ufc_code["create_sub_element"] = "return new %s();" % ufc_code["classname"]
    else:
        cases = ["return new %s_%d();" % (ufc_code["classname"], i) for i in range(num_sub_elements)]
        ufc_code["create_sub_element"] = _generate_switch("i", cases, "return 0;")

    if code_section == "combined":
        return _generate_code(finite_element_combined, ufc_code, options, format)
    elif code_section == "header":
        return _generate_code(finite_element_header, ufc_code, options, format)
    elif code_section == "implementation":
        return _generate_code(finite_element_implementation, ufc_code, options, format)

def _generate_dof_map(code, form_data, options, prefix, i, label, code_section, format):
    "Generate code for ufc::dof_map"

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = format["classname dof_map"](prefix, i, label)

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
    ufc_code["needs_mesh_entities"] = _generate_switch("d", cases, "return false;")

    # Generate code for init_mesh
    ufc_code["init_mesh"] = "__global_dimension = %s;\nreturn false;" % code["global_dimension"]

    # Generate code for init_cell
    ufc_code["init_cell"] = "// Do nothing"

    # Generate code for init_cell_finalize
    ufc_code["init_cell_finalize"] = "// Do nothing"

    # Generate code for global_dimension
    ufc_code["global_dimension"] = "return __global_dimension;"

    # Generate code for local_dimension
    ufc_code["local_dimension"] = "return %s;" % code["local_dimension"]

    # Generate code for max_local_dimension
    ufc_code["max_local_dimension"] = "return %s;" % code["local_dimension"]

    # Generate code for geometric_dimension
    ufc_code["geometric_dimension"] = "return %s;" % code["geometric_dimension"]

    # Generate code for num_facet_dofs
    ufc_code["num_facet_dofs"] = "return %s;" % code["num_facet_dofs"]

    # Generate code for num_entity_dofs (introduced in UFC 1.1)
    ufc_code["num_entity_dofs"] = format["exception"]("Not implemented (introduced in UFC v1.1).")

    # Generate code for tabulate_dofs
    ufc_code["tabulate_dofs"] = _generate_body(code["tabulate_dofs"])

    # Generate code for tabulate_facet_dofs
    ufc_code["tabulate_facet_dofs"] = _generate_switch("facet", [_generate_body(case) for case in code["tabulate_facet_dofs"]])

    # Generate code for tabulate_entity_dofs (introduced in UFC 1.1)
    ufc_code["tabulate_entity_dofs"] = format["exception"]("Not implemented (introduced in UFC v1.1).")

    # Generate code for tabulate_coordinates
    ufc_code["tabulate_coordinates"] = _generate_body(code["tabulate_coordinates"])

    # Generate code for num_sub_dof_maps
    ufc_code["num_sub_dof_maps"] = "return %s;" % code["num_sub_dof_maps"]

    # Generate code for create_sub_dof_map
    num_sub_dof_maps = eval(code["num_sub_dof_maps"])
    if num_sub_dof_maps == 1:
        ufc_code["create_sub_dof_map"] = "return new %s();" % ufc_code["classname"]
    else:
        cases = ["return new %s_%d();" % (ufc_code["classname"], i) for i in range(num_sub_dof_maps)]
        ufc_code["create_sub_dof_map"] = _generate_switch("i", cases, "return 0;")

    if code_section == "combined":
        return _generate_code(dof_map_combined, ufc_code, options, format)
    elif code_section == "header":
        return _generate_code(dof_map_header, ufc_code, options, format)
    elif code_section == "implementation":
        return _generate_code(dof_map_implementation, ufc_code, options, format)

def _generate_form(code, form_data, options, prefix, i, code_section, format):
    "Generate code for ufc::form"

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = format["classname form"](prefix, i)

    # Generate code for members
    ufc_code["members"] = ""

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for signature
    ufc_code["signature"] = "return \"%s\";" % _generate_body(code["signature"])

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
    num_cases = form_data.rank + form_data.num_coefficients
    cases = ["return new %s();" % format["classname finite_element"](prefix, i, (j,)) for j in range(num_cases)]
    ufc_code["create_finite_element"] = _generate_switch("i", cases, "return 0;")

    # Generate code for create_dof_map
    num_cases = form_data.rank + form_data.num_coefficients
    cases = ["return new %s();" % format["classname dof_map"](prefix, i, (j,)) for j in range(num_cases)]
    ufc_code["create_dof_map"] = _generate_switch("i", cases, "return 0;")

    # Generate code for cell_integral
    num_cases = form_data.num_cell_domains
    cases = ["return new %s();" % format["classname cell_integral"](prefix, i, j) for j in range(num_cases)]
    ufc_code["create_cell_integral"] = _generate_switch("i", cases, "return 0;")

    # Generate code for exterior_facet_integral
    num_cases = form_data.num_exterior_facet_domains
    cases = ["return new %s();" % format["classname exterior_facet_integral"](prefix, i, j) for j in range(num_cases)]
    ufc_code["create_exterior_facet_integral"] = _generate_switch("i", cases, "return 0;")

    # Generate code for interior_facet_integral
    num_cases = form_data.num_interior_facet_domains
    cases = ["return new %s();" % format["classname interior_facet_integral"](prefix, i, j) for j in range(num_cases)]
    ufc_code["create_interior_facet_integral"] = _generate_switch("i", cases, "return 0;")

    if code_section == "combined":
        return _generate_code(form_combined, ufc_code, options, format)
    elif code_section == "header":
        return _generate_code(form_header, ufc_code, options, format)
    elif code_section == "implementation":
        return _generate_code(form_implementation, ufc_code, options, format)

def _generate_dolfin_wrappers(generated_forms, prefix, options, format):
    "Generate code for DOLFIN wrappers"

    info("Writing DOLFIN wrappers.")

    # Special case: single element
    if len(generated_forms) == 1 and isinstance(generated_forms[0][1], ElementData):
        fn = UFCFormNames("0",
                          [],
                          format["classname form"](prefix, 0),
                          [format["classname finite_element"](prefix, 0, (0,))],
                          [format["classname dof_map"](prefix, 0, (0,))])
        return generate_dolfin_code(prefix, "", [fn], (0, 0), False) + "\n\n"

    # Generate name data for each form
    form_names = []
    for (i, (form_code, form_data)) in enumerate(generated_forms):
        n = form_data.rank + form_data.num_coefficients
        fn = UFCFormNames("%d" % i,
                          [c.name for c in form_data.coefficients],
                          format["classname form"](prefix, i),
                          [format["classname finite_element"](prefix, i, (j,)) for j in range(n)],
                          [format["classname dof_map"](prefix, i, (j,)) for j in range(n)])
        form_names.append(fn)

    # Collect element signatures
    element_signatures = []
    for (i, (form_code, form_data)) in enumerate(generated_forms):
        element_signatures += [element.__repr__() for element in form_data.elements[:form_data.rank]]

    # Extract common element if any
    if len(element_signatures) > 0 and element_signatures[1:] == element_signatures[:-1]:
        common_space = (0, 0)
    else:
        common_space = None

    # Call DOLFIN wrapper generator
    return generate_dolfin_code(prefix, "", form_names, common_space, False) + "\n\n"

def _generate_jacobian(cell_dimension, integral_type):
    "Generate code for computing jacobian"

    # Choose space dimension
    if cell_dimension == 1:
        jacobian = jacobian_1D
        facet_determinant = facet_determinant_1D
    elif cell_dimension == 2:
        jacobian = jacobian_2D
        facet_determinant = facet_determinant_2D
    else:
        jacobian = jacobian_3D
        facet_determinant = facet_determinant_3D

    # Check if we need to compute more than one Jacobian
    if integral_type == "cell":
        code  = jacobian % {"restriction":  ""}
        code += "\n\n"
        code += scale_factor
    elif integral_type == "exterior facet":
        code  = jacobian % {"restriction":  ""}
        code += "\n\n"
        code += facet_determinant % {"restriction": "", "facet" : "facet"}
    elif integral_type == "interior facet":
        code  = jacobian % {"restriction": choose_map["+"]}
        code += "\n\n"
        code += jacobian % {"restriction": choose_map["-"]}
        code += "\n\n"
        code += facet_determinant % {"restriction": choose_map["+"], "facet": "facet0"}

    return code

def _generate_normal(cell_dimension, integral_type, reference_normal=False):
    "Generate code for computing normal"

    # Choose space dimension
    if cell_dimension == 1:
        normal_direction = normal_direction_1D
        facet_normal = facet_normal_1D
    elif cell_dimension == 2:
        normal_direction = normal_direction_2D
        facet_normal = facet_normal_2D
    else:
        normal_direction = normal_direction_3D
        facet_normal = facet_normal_3D

    if integral_type == "exterior facet":
        code = normal_direction % {"restriction": "", "facet" : "facet"}
        code += facet_normal % {"direction" : "", "restriction": ""}
    elif integral_type == "interior facet":
        code = normal_direction % {"restriction": choose_map["+"], "facet": "facet0"}
        code += facet_normal % {"direction" : "", "restriction": choose_map["+"]}
        code += facet_normal % {"direction" : "!", "restriction": choose_map["-"]}
    return code

def _generate_switch(variable, cases, default=""):
    "Generate switch statement from given variable and cases"

    # Special case: no cases
    if len(cases) == 0:
        return default

    # Special case: one case
    if len(cases) == 1:
        return cases[0]

    # Create switch
    code = "switch (%s)\n{\n" % variable
    for i in range(len(cases)):
        code += "case %d:\n%s\n  break;\n" % (i, indent(cases[i], 2))
    code += "}"
    if not default == "":
        code += "\n" + default

    return code

def _generate_cell_integral(code, form_data, options, prefix, i, label, code_section, format):
    "Generate code for ufc::cell_integral."

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = format["classname cell_integral"](prefix, i, label)

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for members
    # Note special handling of <form prefix> not known at code generation stage!
    body = _generate_body(code["members"])
    body = body.replace("<form prefix>", format["form prefix"](prefix, i))
    ufc_code["members"] = body

    # Generate code for tabulate_tensor
    ufc_code["tabulate_tensor"] = _generate_body(code["tabulate_tensor"])

    if code_section == "combined":
        return _generate_code(cell_integral_combined, ufc_code, options, format)
    elif code_section == "header":
        return _generate_code(cell_integral_header, ufc_code, options, format)
    elif code_section == "implementation":
        return _generate_code(cell_integral_implementation, ufc_code, options, format)

def _generate_exterior_facet_integral(code, form_data, options, prefix, i, label, code_section, format):
    "Generate code for ufc::exterior_facet_integral."

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = format["classname exterior_facet_integral"](prefix, i, label)

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for members
    # Note special handling of <form prefix> not known at code generation stage!
    body = _generate_body(code["members"])
    body = body.replace("<form prefix>", format["form prefix"](prefix, i))
    ufc_code["members"] = body

    # Generate code for tabulate_tensor
    if isinstance(code["tabulate_tensor"], tuple):
        body_lines, cases = code["tabulate_tensor"]
        body = _generate_body(body_lines)
        if len(cases) > 0:
            body += "\n"
            body += _generate_switch("facet", [_generate_body(case) for case in cases])
        ufc_code["tabulate_tensor"] = body
    else:
        ufc_code["tabulate_tensor"] = _generate_body(code["tabulate_tensor"])

    if code_section == "combined":
        return _generate_code(exterior_facet_integral_combined, ufc_code, options, format)
    elif code_section == "header":
        return _generate_code(exterior_facet_integral_header, ufc_code, options, format)
    elif code_section == "implementation":
        return _generate_code(exterior_facet_integral_implementation, ufc_code, options, format)

def _generate_interior_facet_integral(code, form_data, options, prefix, i, label, code_section, format):
    "Generate code for ufc::interior_facet_integral."

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = format["classname interior_facet_integral"](prefix, i, label)

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for members
    # Note special handling of <form prefix> not known at code generation stage!
    body = _generate_body(code["members"])
    body = body.replace("<form prefix>", format["form prefix"](prefix, i))
    ufc_code["members"] = body

    # Generate code for tabulate_tensor
    if isinstance(code["tabulate_tensor"], tuple):
        body_lines, cases = code["tabulate_tensor"]
        body = _generate_body(body_lines)
        if len(cases) > 0:
            body += "\n"
            body += _generate_switch("facet", [_generate_body(case) for case in cases])
        ufc_code["tabulate_tensor"] = body
    else:
        ufc_code["tabulate_tensor"] = _generate_body(code["tabulate_tensor"])

    # Generate code for tabulate_tensor
    if isinstance(code["tabulate_tensor"], tuple):
        body_lines, cases = code["tabulate_tensor"]
        body = _generate_body(body_lines)
        if len(cases) > 0:
            body += "\n"
            body += _generate_switch("facet0", [_generate_switch("facet1", [_generate_body(c) for c in case]) for case in cases])
        ufc_code["tabulate_tensor"] = body
    else:
        ufc_code["tabulate_tensor"] = _generate_body(code["tabulate_tensor"])

    if code_section == "combined":
        return _generate_code(interior_facet_integral_combined, ufc_code, options, format)
    elif code_section == "header":
        return _generate_code(interior_facet_integral_header, ufc_code, options, format)
    elif code_section == "implementation":
        return _generate_code(interior_facet_integral_implementation, ufc_code, options, format)

def _generate_code(format_string, code, options, format):
    "Generate code according to format string and code dictionary"

    # Fix indentation
    for key in code:
        flag = "no-" + key
        if flag in options and options[flag]:
            code[key] = format["exception"]("// Function %s not generated (compiled with -fno-%s)" % (key, key))
        if not key in ["classname", "members"]:
            code[key] = indent(code[key], 4)

    # Generate code
    return format_string % code
