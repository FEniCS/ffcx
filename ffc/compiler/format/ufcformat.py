"Code generation for the UFC 1.0 format"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-08 -- 2009-03-08"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard 2008
# Modified by Dag Lindbo, 2008
# Modified by Johan hake, 2009
# Python modules
import os

# UFC code templates
from ufc_utils import *

# FFC common modules
from ffc.common.utils import *
from ffc.common.log import debug
from ffc.common.constants import *

# FFC format modules
from codesnippets import *
from removeunused import *
from dolfintemplates import *

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
            raise RuntimeError, "Don't know how to compile code for format '%s'." % output_format

        # Attach format
        self.format = {
# operators
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
           "table declaration": "const static double ",
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
           "element tensor": lambda i, k: "A[%d]" % k, # FIXME: Remove argument i, not used
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
           "coordinate map": lambda d: eval("map_coordinates_%dD" % d),
           "facet sign": lambda e: "sign_facet%d" % e,
           "snippet facet signs": lambda d: eval("facet_sign_snippet_%dD" % d),
           "snippet dof map": evaluate_basis_dof_map,
           "snippet eta_interval": eta_interval_snippet,
           "snippet eta_triangle": eta_triangle_snippet,
           "snippet eta_tetrahedron": eta_tetrahedron_snippet,
           "snippet jacobian": lambda d: eval("jacobian_%dD" % d),
           "snippet only jacobian": lambda d: eval("only_jacobian_%dD" % d),

           "snippet combinations": combinations_snippet,
           "snippet transform": lambda d: eval("transform%dD_snippet" % d),
#           "snippet inverse 2D": inverse_jacobian_2D,
#           "snippet inverse 3D": inverse_jacobian_3D,
           "snippet evaluate_dof": lambda d : eval("evaluate_dof_%dD" % d),
           "snippet map_onto_physical": lambda d : eval("map_onto_physical_%dD" % d),
#           "snippet declare_representation": declare_representation,
#           "snippet delete_representation": delete_representation,
           "snippet calculate dof": calculate_dof,
           "get cell vertices" : "const double * const * x = c.coordinates;",
           "generate jacobian": lambda d, i: self.__generate_jacobian(d, i),
           "generate body": lambda d: _generate_body(d),
# misc
           "comment": lambda v: "// %s" % v,
           "pointer": "*",
           "new": "new ",
           "delete": "delete ",
           "cell shape": lambda i: {1: "ufc::interval", 2: "ufc::triangle", 3: "ufc::tetrahedron"}[i],
           "psi index names": {0: lambda i: "f%s" %(i), 1: lambda i: "p%s" %(i),\
                               2: lambda i: "s%s" %(i), 4: lambda i: "fu%s" %(i),\
                               5: lambda i: "pj%s" %(i), 6: lambda i: "c%s" %(i),\
                               7: lambda i: "a%s" %(i)}}

        # Set number of digits for floating point and machine precision
        precision = int(options["precision"])
        f1 = "%%.%dg" % precision
        f2 = "%%.%de" % precision
        def floating_point(v):
            if abs(v) < 100.0:
                return f1 % v
            else:
                return f2 % v
        self.format["floating point"] = floating_point
        self.format["epsilon"] = 10.0*eval("1e-%s" % precision)

    def write(self, generated_forms, prefix, options):
        "Generate UFC 1.0 code for a given list of pregenerated forms"
        debug("Generating code for UFC 1.0")
          
        # Strip directory names from prefix and add output directory
        prefix = prefix.split(os.path.join(' ',' ').split()[0])[-1]
        full_prefix = os.path.join(options["output_dir"], prefix)

        # Generate code for header
        output = ""
        if self.output_format == "ufc":
            output += self.generate_header(prefix, options)
        elif self.output_format == "dolfin":
            output += self.generate_dolfin_header(prefix, options)
        output += "\n"

        if not options["split_implementation"]:

            if self.output_format == "dolfin":
                debug("Generating code for UFC 1.0 with DOLFIN wrappers")

                # Generate UFC code
                output += self.generate_ufc(generated_forms, "UFC_" + prefix, options, "combined")

                # Generate code for DOLFIN wrappers
                output += self._generate_dolfin_wrappers(generated_forms, prefix, options)

            elif self.output_format == "ufc":
                # Generate UFC code
                output += self.generate_ufc(generated_forms, prefix, options, "combined")
            
            # Generate code for footer
            output += self.generate_footer(prefix, options)

            # Write file
            filename = "%s.h" % full_prefix
            file = open(filename, "w")
            file.write(output)
            file.close()
            debug("Output written to " + filename)

        else:

            if self.output_format == "dolfin":
                debug("Generating code for UFC 1.0 with DOLFIN wrappers, split header and implementation")

                # Generate UFC header code
                output += self.generate_ufc(generated_forms, "UFC_" + prefix, options, "header")

                # Generate code for DOLFIN wrappers
                output += self._generate_dolfin_wrappers(generated_forms, prefix, options)

            elif self.output_format == "ufc":
                # Generate UFC code
                output += self.generate_ufc(generated_forms, prefix, options, "header")

            # Generate code for footer
            output += self.generate_footer(prefix, options)

            # Write file
            filename = "%s.h" % full_prefix
            file = open(filename, "w")
            file.write(output)
            file.close()
            debug("Output written to " + filename)

            output = ""
            
            # Generate UFC implementation code
            output += "#include \"%s.h\"\n" % prefix

            if self.output_format == "dolfin":
                output += self.generate_ufc(generated_forms, "UFC_" + prefix, options, "implementation")

            elif self.output_format == "ufc":
                output += self.generate_ufc(generated_forms, prefix, options, "implementation")

            # Write file
            filename = "%s.cpp" % full_prefix
            file = open(filename, "w")
            file.write(output)
            file.close()
            debug("Output written to " + filename)

    def generate_header(self, prefix, options):
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

    def generate_dolfin_header(self, prefix, options):
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

    def generate_footer(self, prefix, options):
        "Generate code for footer"
        return """\
#endif
"""

    def generate_ufc(self, generated_forms, prefix, options, code_section):
        "Generate code for body"

        output = ""
        
        # Iterate over forms
        for i in range(len(generated_forms)):

            # Get pregenerated code, form data and prefix
            (form_code, form_data) = generated_forms[i]
            form_prefix = self.compute_prefix(prefix, generated_forms, i, options)

            # Generate code for ufc::finite_element(s)
            for (label, sub_element) in form_code["finite_elements"]:
                output += self.__generate_finite_element(sub_element, form_data, options, form_prefix, label, code_section)
                output += "\n"

            # Generate code for ufc::dof_map(s)
            for (label, sub_dof_map) in form_code["dof_maps"]:
                output += self.__generate_dof_map(sub_dof_map, form_data, options, form_prefix, label, code_section)
                output += "\n"

            # Generate code for ufc::cell_integral
            if form_code.has_key("cell_integrals"):
                for (label, code) in form_code["cell_integrals"]:
                    output += _generate_cell_integral(code, form_data, options, form_prefix, label, code_section, self.format["exception"])
                    output += "\n"

            # Generate code for ufc::exterior_facet_integral
            if form_code.has_key("exterior_facet_integrals"):
                for (label, code) in form_code["exterior_facet_integrals"]:
                    output += _generate_exterior_facet_integral(code, form_data, options, form_prefix, label, code_section, self.format["exception"])
                    output += "\n"

            # Generate code for ufc::interior_facet_integral
            if form_code.has_key("interior_facet_integrals"):
                for (label, code) in form_code["interior_facet_integrals"]:
                    output += _generate_interior_facet_integral(code, form_data, options, form_prefix, label, code_section, self.format["exception"])
                    output += "\n"

            # Generate code for ufc::form
            if "form" in form_code:
                output += self.__generate_form(form_code["form"], form_data, options, form_prefix, code_section)
                output += "\n"

        return output

    def compute_prefix(self, prefix, generated_forms, i, options):
        "Compute prefix for form i"

        # Get form ranks
        ranks = [form_data.rank for (form_code, form_data) in generated_forms]

        # Return prefixFunctional, prefixLinearForm or prefixBilinearForm
        # when we have exactly one form of ranks 0, 1 or 2
        count = [ranks.count(0), ranks.count(1), ranks.count(2)]
        if len(ranks) <= 3 and sum(count) > 0 and min(count) >= 0 and max(count) <= 1 and options["form_postfix"]:
            postfixes = ["Functional", "LinearForm", "BilinearForm"]
            return "%s%s" % (prefix, postfixes[ranks[i]])

        # Return prefix_i if we have more than one rank
        if len(ranks) > 1:
            return "%s_%d" % (prefix, i)

        # Else, just return prefix
        return prefix

    def __generate_finite_element(self, code, form_data, options, prefix, label, code_section):
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
        ufc_code["value_dimension"] = self.__generate_switch("i", cases, "return 0;")

        # Generate code for evaluate_basis (and vectorised counterpart)
        ufc_code["evaluate_basis"] = _generate_body(code["evaluate_basis"])
        ufc_code["evaluate_basis_all"] = _generate_body(code["evaluate_basis_all"])

        # Generate code for evaluate_basis_derivatives (and vectorised counterpart)
        ufc_code["evaluate_basis_derivatives"] = _generate_body(code["evaluate_basis_derivatives"])
        ufc_code["evaluate_basis_derivatives_all"] = _generate_body(code["evaluate_basis_derivatives_all"])

        # Generate code for evaluate_dof
        ufc_code["evaluate_dof"] = _generate_body(code["evaluate_dof"])

        # Generate code for evaluate_dofs (introduced in UFC 1.1)
        ufc_code["evaluate_dofs"] = self.format["exception"]("Not implemented (introduced in UFC v1.1).")

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
            ufc_code["create_sub_element"] = self.__generate_switch("i", cases, "return 0;")
        
        if code_section == "combined":
            return _generate_code(finite_element_combined, ufc_code, options, self.format["exception"])
        elif code_section == "header":
            return _generate_code(finite_element_header, ufc_code, options, self.format["exception"])
        elif code_section == "implementation":
            return _generate_code(finite_element_implementation, ufc_code, options, self.format["exception"])

    def __generate_dof_map(self, code, form_data, options, prefix, label, code_section):
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
        ufc_code["needs_mesh_entities"] = self.__generate_switch("d", cases, "return false;")

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

        # Generate code for geometric_dimension
        ufc_code["geometric_dimension"] = "return %s;" % code["geometric_dimension"]

        # Generate code for num_facet_dofs
        ufc_code["num_facet_dofs"] = "return %s;" % code["num_facet_dofs"]

        # Generate code for num_entity_dofs (introduced in UFC 1.1)
        ufc_code["num_entity_dofs"] = self.format["exception"]("Not implemented (introduced in UFC v1.1).")

        # Generate code for tabulate_dofs
        ufc_code["tabulate_dofs"] = _generate_body(code["tabulate_dofs"])

        # Generate code for tabulate_facet_dofs
        ufc_code["tabulate_facet_dofs"] = self.__generate_switch("facet", [_generate_body(case) for case in code["tabulate_facet_dofs"]])

        # Generate code for tabulate_entity_dofs (introduced in UFC 1.1)
        ufc_code["tabulate_entity_dofs"] = self.format["exception"]("Not implemented (introduced in UFC v1.1).")

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
            ufc_code["create_sub_dof_map"] = self.__generate_switch("i", cases, "return 0;")

        if code_section == "combined":
            return _generate_code(dof_map_combined, ufc_code, options, self.format["exception"])
        elif code_section == "header":
            return _generate_code(dof_map_header, ufc_code, options, self.format["exception"])
        elif code_section == "implementation":
            return _generate_code(dof_map_implementation, ufc_code, options, self.format["exception"])

    def __generate_form(self, code, form_data, options, prefix, code_section):
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
        cases = ["return new %s_finite_element_%d();" % (prefix, i) for i in range(num_cases)]
        ufc_code["create_finite_element"] = self.__generate_switch("i", cases, "return 0;")

        # Generate code for create_dof_map
        num_cases = form_data.rank + form_data.num_coefficients
        cases = ["return new %s_dof_map_%d();" % (prefix, i) for i in range(num_cases)]
        ufc_code["create_dof_map"] = self.__generate_switch("i", cases, "return 0;")

        # Generate code for cell_integral
        num_cases = form_data.num_cell_integrals
        cases = ["return new %s_cell_integral_%d();" % (prefix, i) for i in range(num_cases)]
        ufc_code["create_cell_integral"] = self.__generate_switch("i", cases, "return 0;")

        # Generate code for exterior_facet_integral
        num_cases = form_data.num_exterior_facet_integrals
        cases = ["return new %s_exterior_facet_integral_%d();" % (prefix, i) for i in range(num_cases)]
        ufc_code["create_exterior_facet_integral"] = self.__generate_switch("i", cases, "return 0;")

        # Generate code for interior_facet_integral
        num_cases = form_data.num_interior_facet_integrals
        cases = ["return new %s_interior_facet_integral_%d();" % (prefix, i) for i in range(num_cases)]
        ufc_code["create_interior_facet_integral"] = self.__generate_switch("i", cases, "return 0;")

        if code_section == "combined":
            return _generate_code(form_combined, ufc_code, options, self.format["exception"])
        elif code_section == "header":
            return _generate_code(form_header, ufc_code, options, self.format["exception"])
        elif code_section == "implementation":
            return _generate_code(form_implementation, ufc_code, options, self.format["exception"])

    def _generate_dolfin_wrappers(self, generated_forms, prefix, options):
        "Generate code for DOLFIN wrappers"

        # FIXME: Use module from dolfin_utils

        # We generate two versions of all constructors, one using references (_r)
        # and one using shared pointers (_s)

        output = dolfin_includes

        # Extract common test space if any
        test_element = None
        elements_string_map = {}
        test_elements = []
        for (form_code, form_data) in generated_forms:
            if form_data.rank >= 1:
                test_elements.append(form_data.ffc_elements[0].__repr__())
                elements_string_map[form_data.ffc_elements[0].__repr__()] = form_data.ffc_elements[0]
        if len(test_elements) > 0 and test_elements[1:] == test_elements[:-1]:
            test_element = test_elements[0]
        elif len(test_elements) > 0:
            raise RuntimeError, "Unable to extract test space (not uniquely defined)."

        # Extract common trial element if any
        trial_element = None
        trial_elements = []
#        trial_elements = [form_data.elements[0] for (form_code, form_data) in generated_forms if form_data.rank >= 1]
        for (form_code, form_data) in generated_forms:
            if form_data.rank >= 1:
                trial_elements.append(form_data.ffc_elements[0].__repr__())
                elements_string_map[form_data.ffc_elements[0].__repr__()] = form_data.ffc_elements[0]
        if len(trial_elements) > 0 and trial_elements[1:] == trial_elements[:-1]:
            trial_element = trial_elements[0]
        elif len(trial_elements) > 0:
            raise RuntimeError, "Unable to extract trial space (not uniquely defined)."

        # Extract common coefficient element if any
        coefficient_element = None
        coefficients = [c for (form_code, form_data) in generated_forms for c in form_data.coefficients]
        coefficient_elements = []
        for c in coefficients:
#            coefficient_elements.append(c.e0)
            coefficient_elements.append(c.element.__repr__())
            elements_string_map[c.element.__repr__()] = c.element
        if len(coefficient_elements) > 0 and coefficient_elements[1:] == coefficient_elements[:-1]:
            coefficient_element = coefficient_elements[0]

        # Extract common element if any
        common_element = None
        elements = test_elements + trial_elements # + coefficient_elements
        if len(elements) > 0 and elements[1:] == elements[:-1]:
            common_element = elements[-1]

        # Build map from elements to forms
        element_map = {}    
        for i in range(len(generated_forms)):
            (form_code, form_data) = generated_forms[i]
            form_prefix = self.compute_prefix(prefix, generated_forms, i, options)
            for j in range(len(form_data.ffc_elements)):
#                element_map[form_data.elements[j]] = (form_prefix, j)
                element_map[form_data.ffc_elements[j].__repr__()] = (form_prefix, j)

        # Generate code for function spaces
        for i in range(len(generated_forms)):
            (form_code, form_data) = generated_forms[i]
            form_prefix = self.compute_prefix(prefix, generated_forms, i, options)
            for j in range(form_data.rank):
#                output += self._generate_function_space(form_data.elements[j],
                output += self._generate_function_space(form_data.ffc_elements[j].__repr__(),
                                                   "%sFunctionSpace%d" % (form_prefix, j),
                                                   element_map)
                output += "\n"
            for j in range(form_data.num_coefficients):
#                output += self._generate_function_space(form_data.elements[form_data.rank + j],
                output += self._generate_function_space(form_data.ffc_elements[form_data.rank + j].__repr__(),
                                                   "%sCoefficientSpace%d" % (form_prefix, j),
                                                   element_map)
                output += "\n"
            # If we compiled just one element
            if form_data.rank < 0:
                if not len(form_data.ffc_elements) == 1:
                    raise RuntimeError("Expected just one element")
                output += self._generate_function_space(form_data.ffc_elements[0].__repr__(),
                                                       "%sFunctionSpace" % (form_prefix),
                                                       element_map)
                output += "\n"

        # Generate code for special function spaces
        if not test_element is None:
            output += self._generate_function_space(test_element, prefix + "TestSpace", element_map) + "\n"
        if not trial_element is None:
            output += self._generate_function_space(trial_element, prefix + "TrialSpace", element_map) + "\n"
        if not coefficient_element is None:
            output += self._generate_function_space(coefficient_element, prefix + "CoefficientSpace", element_map) + "\n"
        if not common_element is None:
            output += self._generate_function_space(common_element, prefix + "FunctionSpace", element_map) + "\n"

        # Generate wrappers for all forms
        element_map = {}
        for i in range(len(generated_forms)):

            (form_code, form_data) = generated_forms[i]
            form_prefix = self.compute_prefix(prefix, generated_forms, i, options)

            # Generate code for coefficient member variables
            coefficient_names = [c.name for c in form_data.coefficients]
            n = len(coefficient_names)
            coefficient_classes = ["%sCoefficient%d" % (form_prefix, j) for j in range(n)]
            if n == 0:
                coefficient_members = ""
            else:
                coefficient_members = "\n  // Coefficients\n" + "\n".join(["  %s %s;" % (coefficient_classes[j], coefficient_names[j]) for j in range(n)]) + "\n"

            # Generate code for initialization of coefficients
            if n == 0:
                coefficient_init = ""
            else:
                coefficient_init = ", " + ", ".join(["%s(*this)" % coefficient_names[j] for j in range(n)])

            # Generate code for coefficient classes
            for j in range(n):
                output += coefficient_class % (coefficient_classes[j],
                                               coefficient_classes[j],
                                               coefficient_classes[j],
                                               coefficient_classes[j],
                                               form_prefix, j, j, coefficient_names[j])

            # Generate constructors
            assign_coefficients = "\n".join(["    this->%s = w%d;" % (coefficient_names[k], k) for k in range(form_data.num_coefficients)])
            if form_data.num_coefficients > 0:
                assign_coefficients += "\n\n"
            constructor_args_r  = ", ".join(["const dolfin::FunctionSpace& V%d" % k for k in range(form_data.rank)])
            constructor_args_rc = ", ".join(["const dolfin::FunctionSpace& V%d" % k for k in range(form_data.rank)] +
                                            ["dolfin::Function& w%d" % k for k in range(form_data.num_coefficients)])
            constructor_args_s  = ", ".join(["boost::shared_ptr<const dolfin::FunctionSpace> V%d" % k for k in range(form_data.rank)])
            constructor_args_sc  = ", ".join(["boost::shared_ptr<const dolfin::FunctionSpace> V%d" % k for k in range(form_data.rank)] +
                                             ["dolfin::Function& w%d" % k for k in range(form_data.num_coefficients)])
            constructor_body_r  = "\n".join([add_function_space_r % (k, k, k) for k in range(form_data.rank)])
            constructor_body_rc = "\n".join([add_function_space_r % (k, k, k) for k in range(form_data.rank)])
            constructor_body_s  = "\n".join([add_function_space_s % k for k in range(form_data.rank)])
            constructor_body_sc = "\n".join([add_function_space_s % k for k in range(form_data.rank)])
            if form_data.rank > 0:
                constructor_body_r  += "\n\n"
                constructor_body_rc += "\n\n"
                constructor_body_s  += "\n\n"
                constructor_body_sc += "\n\n"
            constructor_body_r  += "\n".join([add_coefficient_r for k in range(form_data.num_coefficients)])
            constructor_body_rc += "\n".join([add_coefficient_r for k in range(form_data.num_coefficients)])
            constructor_body_s  += "\n".join([add_coefficient_s for k in range(form_data.num_coefficients)])
            constructor_body_sc += "\n".join([add_coefficient_s for k in range(form_data.num_coefficients)])
            if form_data.num_coefficients > 0:
                constructor_body_r  += "\n\n"
                constructor_body_rc += "\n\n"
                constructor_body_s  += "\n\n"
                constructor_body_sc += "\n\n"
            constructor_body_rc += assign_coefficients
            constructor_body_sc += assign_coefficients
            constructor_body_r  += "    _ufc_form = boost::shared_ptr<const ufc::form>(new UFC_%s());" % form_prefix
            constructor_body_rc += "    _ufc_form = boost::shared_ptr<const ufc::form>(new UFC_%s());" % form_prefix
            constructor_body_s  += "    _ufc_form = boost::shared_ptr<const ufc::form>(new UFC_%s());" % form_prefix
            constructor_body_sc += "    _ufc_form = boost::shared_ptr<const ufc::form>(new UFC_%s());" % form_prefix

            # Generate class in different ways depending on the situation
            if form_data.rank > 0:
                if form_data.num_coefficients > 0:
                    output += form_class_vc % (form_prefix,
                                               form_prefix, constructor_args_r,  coefficient_init, constructor_body_r,
                                               form_prefix, constructor_args_s,  coefficient_init, constructor_body_s,
                                               form_prefix, constructor_args_rc, coefficient_init, constructor_body_rc,
                                               form_prefix, constructor_args_sc, coefficient_init, constructor_body_sc,
                                               form_prefix, coefficient_members)
                else:
                    output += form_class_v % (form_prefix,
                                              form_prefix, constructor_args_r,  coefficient_init, constructor_body_r,
                                              form_prefix, constructor_args_s,  coefficient_init, constructor_body_s,
                                              form_prefix, coefficient_members)
            else:
                if form_data.num_coefficients > 0:
                    output += form_class_c % (form_prefix,
                                              form_prefix, constructor_args_r,  coefficient_init, constructor_body_r,
                                              form_prefix, constructor_args_rc, coefficient_init, constructor_body_rc,
                                              form_prefix, coefficient_members)
#                else:
#                    output += form_class % (form_prefix,
#                                            form_prefix, constructor_args_r,  coefficient_init, constructor_body_r,
#                                            form_prefix, coefficient_members)

        return output

    def _generate_function_space(self, element, classname, element_map):
        "Generate code for function space"

        function_space_class = """\
class %s : public dolfin::FunctionSpace
{
public:

  %s(const dolfin::Mesh& mesh)
    : dolfin::FunctionSpace(boost::shared_ptr<const dolfin::Mesh>(&mesh, dolfin::NoDeleter<const dolfin::Mesh>()),
                            boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new %s()))),
                            boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new %s()), mesh)))
  {
    // Do nothing
  }

};
"""
        (form_prefix, element_number) = element_map[element]
        element_class = "UFC_%s_finite_element_%d" % (form_prefix, element_number)
        dofmap_class = "UFC_%s_dof_map_%d" % (form_prefix, element_number)
        return function_space_class % (classname, classname, element_class, dofmap_class)

    def __generate_jacobian(self, cell_dimension, integral_type):
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

    def __generate_switch(self, variable, cases, default=""):
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

def _generate_cell_integral(code, form_data, options, prefix, postfix, code_section, format_exception):
    "Generate code for ufc::cell_integral."

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = "%s_cell_integral_%s" % (prefix, postfix)

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for members
    # Note special handling of <form prefix> not known at code generation stage!
    body = _generate_body(code["members"])
    body = body.replace("<form prefix>", prefix)
    ufc_code["members"] = body

    # Generate code for tabulate_tensor
    ufc_code["tabulate_tensor"] = _generate_body(code["tabulate_tensor"])

    if code_section == "combined":
        return _generate_code(cell_integral_combined, ufc_code, options, format_exception)
    elif code_section == "header":
        return _generate_code(cell_integral_header, ufc_code, options, format_exception)
    elif code_section == "implementation":
        return _generate_code(cell_integral_implementation, ufc_code, options, format_exception)

def _generate_exterior_facet_integral(code, form_data, options, prefix, postfix, code_section, format_exception):
    "Generate code for ufc::exterior_facet_integral."

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = "%s_exterior_facet_integral_%s" % (prefix, postfix)

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for members
    # Note special handling of <form prefix> not known at code generation stage!
    body = _generate_body(code["members"])
    body = body.replace("<form prefix>", prefix)
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
        return _generate_code(exterior_facet_integral_combined, ufc_code, options, format_exception)
    elif code_section == "header":
        return _generate_code(exterior_facet_integral_header, ufc_code, options, format_exception)
    elif code_section == "implementation":
        return _generate_code(exterior_facet_integral_implementation, ufc_code, options, format_exception)

def _generate_interior_facet_integral(code, form_data, options, prefix, postfix, code_section, format_exception):
    "Generate code for ufc::interior_facet_integral."

    ufc_code = {}

    # Set class name
    ufc_code["classname"] = "%s_interior_facet_integral_%s" % (prefix, postfix)

    # Generate code for constructor
    ufc_code["constructor"] = "// Do nothing"

    # Generate code for destructor
    ufc_code["destructor"] = "// Do nothing"

    # Generate code for members
    # Note special handling of <form prefix> not known at code generation stage!
    body = _generate_body(code["members"])
    body = body.replace("<form prefix>", prefix)
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
        return _generate_code(interior_facet_integral_combined, ufc_code, options, format_exception)
    elif code_section == "header":
        return _generate_code(interior_facet_integral_header, ufc_code, options, format_exception)
    elif code_section == "implementation":
        return _generate_code(interior_facet_integral_implementation, ufc_code, options, format_exception)

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

def _generate_code(format_string, code, options, format_exception):
    "Generate code according to format string and code dictionary"

    # Fix indentation
    for key in code:
        flag = "no-" + key
        if flag in options and options[flag]:
            code[key] = format_exception("// Function %s not generated (compiled with -fno-%s)" % (key, key))
        if not key in ["classname", "members"]:
            code[key] = indent(code[key], 4)

    # Generate code
    return format_string % code
