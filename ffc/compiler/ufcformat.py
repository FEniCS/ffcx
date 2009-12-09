"Code generation for the UFC 1.0 format"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-08"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.
# Modified by Dag Lindbo, 2008.
# Modified by Johan Hake, 2009.
# Modified by Garth N. Wells, 2009.
# Last changed: 2009-12-09

# Python modules
import os
import platform

# UFC code templates
# UFC build
from ufc_utils import build_ufc_module

# UFC function
from ufc_utils import function_combined
from ufc_utils import function_header
from ufc_utils import function_implementation

# UFC finite_element
from ufc_utils import finite_element_combined
from ufc_utils import finite_element_header
from ufc_utils import finite_element_implementation

# UFC dof_map
from ufc_utils import dof_map_combined
from ufc_utils import dof_map_header
from ufc_utils import dof_map_implementation

# UFC integrals
from ufc_utils import cell_integral_combined
from ufc_utils import cell_integral_header
from ufc_utils import cell_integral_implementation
from ufc_utils import exterior_facet_integral_combined
from ufc_utils import exterior_facet_integral_header
from ufc_utils import exterior_facet_integral_implementation
from ufc_utils import interior_facet_integral_combined
from ufc_utils import interior_facet_integral_header
from ufc_utils import interior_facet_integral_implementation

# UFC form
from ufc_utils import form_combined
from ufc_utils import form_header
from ufc_utils import form_implementation

# DOLFIN wrapper generator
try:
    from dolfin_utils.wrappers import generate_dolfin_code
    from dolfin_utils.wrappers import UFCFormNames
    dolfin_utils_imported = True
except:
    dolfin_utils_imported = False

# FFC common modules
#from ffc.common.utils import *
from ffc.common.log import info, error
from ffc.common.constants import FFC_VERSION

# FFC format modules
# TODO: Finish this import list
#from codesnippets import evaluate_basis_dof_map, eta_interval_snippet, eta_triangle_snippet, eta_tetrahedron_snippet, combinations_snippet, calculate_dof, map_coordinates_interval, map_coordinates_triangle, map_coordinates_tetrahedron, transform_interval_snippet, transform_triangle_snippet, transform_tetrahedron_snippet, map_onto_physical_1D, map_onto_physical_2D, map_onto_physical_3D, jacobian_1D, jacobian_2D, jacobian_3D, facet_determinant_1D, facet_determinant_2D, facet_determinant_3D, scale_factor, normal_direction_2D

from codesnippets import *

# FFC codegeneration common modules
from codeutils import indent
from removeunused import remove_unused
from compiler import ElementData

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
        self.format = {
            #
            # Operators
            #
            "times equal": lambda i, j: "%s *= %s;" %(i,j),
            "add equal": lambda i, j: "%s += %s;" % (i,j),
            "inverse": lambda v: "(1.0/%s)" % v,
            "absolute value": lambda v: "std::abs(%s)" % v,
            "sqrt": lambda v: "std::sqrt(%s)" % v,
            "add": lambda v: " + ".join(v),
            "subtract": lambda v: " - ".join(v),
            "multiply": lambda v: "*".join(v),
            "division": "/",
            "power": lambda base, exp: power_options[exp >= 0](self.format["multiply"]([str(base)]*abs(exp))),
            "std power": lambda base, exp: "std::pow(%s, %s)" % (base, exp),
            "exp": lambda v: "std::exp(%s)" % v,
            "ln": lambda v: "std::log(%s)" % v,
            "cos": lambda v: "std::cos(%s)" % v,
            "sin": lambda v: "std::sin(%s)" % v,
            # bool operators
            "logical and": " && ",
            "logical or": " || ",
            "is equal": " == ",
            "not equal": " != ",
            "less than": " < ",
            "greater than": " > ",
            "bool": lambda v: {True: "true", False: "false"}[v],
            # formating
            "floating point": lambda v: "<not defined>",
            "epsilon": "<not defined>",
            "grouping": lambda v: "(%s)" % v,
            "block": lambda v: "{%s}" % v,
            "block begin": "{",
            "block end": "}",
            "separator": ", ",
            "block separator": ",\n",
            #           "block separator": ",",
            "new line": "\\\n",
            "end line": ";",
            "space": " ",
            # IO
            "exception": lambda v: "throw std::runtime_error(\"%s\");" % v,
            # declarations
            "float declaration": "double ",
            "const float declaration": "const double ",
            "static float declaration": "static double ",
            "uint declaration": "unsigned int ",
            "const uint declaration": "const unsigned int ",
            "static const uint declaration": "static const unsigned int ",
            "static uint declaration": "static unsigned int ",
            "table declaration": "static const double ",
            # variable names
            "element tensor quad": "A",
            "integration points": "ip",
            "first free index": "j",
            "second free index": "k",
            "free secondary indices":["r","s","t","u"],
            "derivatives": lambda i,j,k,l: "dNdx%d_%d[%s][%s]" % (i,j,k,l),
            "element coordinates": lambda i,j: "x[%s][%s]" % (i,j),
            "weight": lambda i: "W%d" % (i),
            "weights": lambda i,j: self.format["weight"](i) + "[%s]" % (j),
            "psis": "P",
            "function value": "F",
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
            "normal component": lambda r, j: "n%s%s" % (choose_map[r], j),
            "tmp declaration": lambda j, k: "const double " + self.format["tmp access"](j, k),
            "tmp access": lambda j, k: "tmp%d_%d" % (j, k),
            "determinant": lambda r: "detJ%s" % choose_map[r],
            "scale factor": "det",
            "constant": lambda j: "c%d" % j,
            "coefficient table": lambda j, k: "w[%d][%d]" % (j, k),
            "coefficient": lambda j, k: "w[%d][%d]" % (j, k),
            "coeff": "w",
            "modified coefficient declaration": lambda i, j, k, l: "const double c%d_%d_%d_%d" % (i, j, k, l),
            "modified coefficient access": lambda i, j, k, l: "c%d_%d_%d_%d" % (i, j, k, l),
            "transform": lambda type, j, k, r: "%s" % (transform_options[type](choose_map[r], j, k)),
            "transform_ufl": lambda type, j, k, r: "%s" % (transform_options_ufl[type](choose_map[r], j, k)),
            "reference tensor" : lambda j, i, a: None,
            "geometry tensor declaration": lambda j, a: "const double " + self.format["geometry tensor access"](j, a),
            "geometry tensor access": lambda j, a: "G%d_%s" % (j, "_".join(["%d" % index for index in a])),
            "geometry tensor": "G",
            "element tensor": lambda i: "A[%d]" % i,
            "sign tensor": lambda type, i, k: "S%s%s_%d" % (type, i, k),
            "sign tensor declaration": lambda s: "const int " + s,
            "signs": "S",
            "vertex values": lambda i: "vertex_values[%d]" % i,
            "dof values": lambda i: "dof_values[%d]" % i,
            "dofs": lambda i: "dofs[%d]" % i,
            "entity index": lambda d, i: "c.entity_indices[%d][%d]" % (d, i),
            "num entities": lambda dim : "m.num_entities[%d]" % dim,
            "offset declaration": "unsigned int offset",
            "offset access": "offset",
            "nonzero columns": lambda i: "nzc%d" % i,
            # access
            "array access": lambda i: "[%s]" %(i),
            "matrix access": lambda i,j: "[%s][%s]" %(i,j),
            "secondary index": lambda i: "_%s" %(i),
            # program flow
            "dof map if": lambda i,j: "if (%d <= %s && %s <= %d)" %(i,\
                                                                    self.format["argument basis num"], self.format["argument basis num"], j),
            "loop": lambda i,j,k: "for (unsigned int %s = %s; %s < %s; %s++)"% (i, j, i, k, i),
            "if": "if",
            "return": lambda v: "return %s;" % v,
            # snippets
            "coordinate map": lambda s: eval("map_coordinates_%s" % s),
            "facet sign": lambda e: "sign_facet%d" % e,
            "snippet facet signs": lambda d: eval("facet_sign_snippet_%dD" % d),
            "snippet dof map": evaluate_basis_dof_map,
            "snippet eta_interval": eta_interval_snippet,
            "snippet eta_triangle": eta_triangle_snippet,
            "snippet eta_tetrahedron": eta_tetrahedron_snippet,
            "snippet jacobian": lambda d: eval("jacobian_%dD" % d),
            "snippet only jacobian": lambda d: eval("only_jacobian_%dD" % d),
            "snippet normal": lambda d: eval("facet_normal_%dD" %d),
            "snippet combinations": combinations_snippet,
            "snippet transform": lambda s: eval("transform_%s_snippet" % s),
            #           "snippet inverse 2D": inverse_jacobian_2D,
            #           "snippet inverse 3D": inverse_jacobian_3D,
            "snippet evaluate_dof": lambda d : eval("evaluate_dof_%dD" % d),
            "snippet map_onto_physical": lambda d : eval("map_onto_physical_%dD" % d),
            #           "snippet declare_representation": declare_representation,
            #           "snippet delete_representation": delete_representation,
            "snippet calculate dof": calculate_dof,
            "get cell vertices" : "const double * const * x = c.coordinates;",
            "generate jacobian": lambda d, i: _generate_jacobian(d, i),
            "generate normal": lambda d, i: _generate_normal(d, i),
            "generate body": lambda d: _generate_body(d),
            # misc
            "comment": lambda v: "// %s" % v,
            "pointer": "*",
            "new": "new ",
            "delete": "delete ",
            "cell shape": lambda i: {"interval": "ufc::interval",
                                     "triangle": "ufc::triangle",
                                     "tetrahedron": "ufc::tetrahedron"}[i],
            "psi index names": {0: lambda i: "f%s" %(i), 1: lambda i: "p%s" %(i),\
                                2: lambda i: "s%s" %(i), 4: lambda i: "fu%s" %(i),\
                                5: lambda i: "pj%s" %(i), 6: lambda i: "c%s" %(i),\
                                7: lambda i: "a%s" %(i)},
            #
            # Class names
            #
            "form prefix": \
            lambda prefix, i: "%s_%d" % (prefix.lower(), i),
            "classname finite_element": \
            lambda prefix, i, label: "%s_%d_finite_element_%s" % (prefix.lower(), i, "_".join([str(j) for j in label])),
            "classname dof_map": \
            lambda prefix, i, label: "%s_%d_dof_map_%s" % (prefix.lower(), i, "_".join([str(j) for j in label])),
            "classname form": \
            lambda prefix, i: "%s_form_%d" % (prefix.lower(), i),
            "classname cell_integral": \
            lambda prefix, i, label: "%s_%d_cell_integral_%s" % (prefix.lower(), i, label),
            "classname interior_facet_integral": \
            lambda prefix, i, label: "%s_%d_interior_facet_integral_%s" % (prefix.lower(), i, label),
            "classname exterior_facet_integral": \
            lambda prefix, i, label: "%s_%d_exterior_facet_integral_%s" % (prefix.lower(), i, label)}

        # Set number of digits for floating point and machine precision
        precision = int(options["precision"])
        f1 = "%%.%dg" % precision
        f2 = "%%.%de" % precision
        def floating_point(v):
            "Format floating point number."
            if abs(v) < 100.0:
                return f1 % v
            else:
                return f2 % v
        def floating_point_windows(v):
            "Format floating point number for Windows (remove extra leading zero in exponents)."
            return floating_point(v).replace("e-0", "e-").replace("e+0", "e+") 
        if platform.system() == "Windows":
            self.format["floating point"] = floating_point_windows
        else:
            self.format["floating point"] = floating_point
        self.format["epsilon"] = 10.0*eval("1e-%s" % precision)

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

def _generate_switch(variable, cases, default = ""):
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

def _generate_body(declarations):
    "Generate function body from list of declarations or statements."

    if not isinstance(declarations, list):
        declarations = [declarations]
    lines = []
    for declaration in declarations:
        if isinstance(declaration, tuple):
            lines += ["%s = %s;" % declaration]
        else:
            lines += ["%s" % declaration]
    return "\n".join(lines)

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

