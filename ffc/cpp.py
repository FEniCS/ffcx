"This module defines rules and algorithms for generating C++ code."

format = {"add": lambda v: _add(v),
          "bool": lambda v: {True: "true", False: "false"}[v],
          "exception": lambda v: "throw std::runtime_error(\"%s\");" % v,
          "return": lambda v: "return %s;" % str(v),
          "multiply": lambda v: _multiply(v),
          "switch": lambda d, cases: _generate_switch(d, cases)}

def _multiply(list):
    non_zero_factors = []
    for e in list:
        if e == "0":
            return ""
        elif e == "1":
            pass
        else:
            non_zero_factors += [e]
    return "*".join(non_zero_factors)

def _add(list):
    return " + ".join([e for e in list if (e != "0" and e != "")])

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

def indent(block, num_spaces):
    "Indent each row of the given string block with n spaces."
    indentation = " " * num_spaces
    return indentation + ("\n" + indentation).join(block.split("\n"))

# FIXME: Major cleanup needed, remove as much as possible
from codesnippets import *

# FIXME: KBO: temporary hack to get dictionary working.
from constants import FFC_OPTIONS
import platform
options=FFC_OPTIONS.copy()

# Old dictionary, move the stuff we need to the new dictionary above
format_old = {
    # Operators
    #
    "times equal": lambda i, j: "%s *= %s;" %(i, j),
    "add equal": lambda i, j: "%s += %s;" % (i, j),
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
    format_old["floating point"] = floating_point_windows
else:
    format_old["floating point"] = floating_point

format_old["epsilon"] = 10.0*eval("1e-%s" % precision)

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
