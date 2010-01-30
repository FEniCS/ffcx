"This module defines rules and algorithms for generating C++ code."

__author__ = "Anders Logg (logg@simula.no) and friends"
__date__ = "2009-12-16"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard 2010
# Modified by Marie E. Rognes 2010
# Last changed: 2010-01-30

# Python modules
import re, numpy, platform

# FFC modules
from ffc.log import debug
from ffc.constants import FFC_OPTIONS

# FIXME: AL: This files needs cleaning up!

# FIXME: AL: In places where we have non-implemented functions
# FIXME: print a warning message instead of throwing an exception
# FIXME: we should throw an exception and instead have a command-line
# FIXME: option for converting exceptions to warnings that can be
# FIXME: used from the regression test script.

# Formatting rules
# FIXME: KBO: format is a builtin_function, i.e., we should use a different name.
format = {}


# Program flow
format.update({"return":      lambda v: "return %s;" % str(v),
               "grouping":    lambda v: "(%s)" % v,
               "block":       lambda v: "{%s}" % v,
               "block begin": "{",
               "block end":   "}",
               "list":        lambda v: format["block"](format["separator"].join([str(l) for l in v])),
               "switch":      lambda v, cases, default=None, numbers=None: _generate_switch(v, cases, default, numbers),
               "exception":   lambda v: "throw std::runtime_error(\"%s\");" % v,
               "warning":     lambda v: 'std::cerr << "*** FFC warning: " << "%s" << std::endl;' % v,
               "comment":     lambda v: "// %s" % v,
               "if":          lambda c, v: "if (%s)\n{\n%s\n}\n" % (c, v),
               "loop":        lambda i, j, k: "for (unsigned int %s = %s; %s < %s; %s++)"% (i, j, i, k, i),
               "is equal": " == ",
               "do nothing":  "// Do nothing"})

# Declarations
format.update({"declaration": lambda t, n, v=None: _declaration(t, n, v),
               "float declaration": "double ",
               "uint declaration": "unsigned int ",
               "static const uint declaration": "static const unsigned int ",
               "static const float declaration": "static const double ",
               "const float declaration":
                   lambda v, w: "const double %s = %s;" % (v, w),
               "const uint declaration":
                   lambda v, w: "const unsigned int %s = %s;" % (v, w)})

# Mathematical operators
format.update({"add":           lambda v: " + ".join(v),
               "iadd":          lambda v, w: "%s += %s;" % (str(v), str(w)),
               "sub":           lambda v: " - ".join(v),
               "mul":           lambda v: "*".join(v),
               "imul":          lambda v, w: "%s *= %s;" % (str(v), str(w)),
               "div":           lambda v, w: "%s/%s" % (str(v), str(w)),
               "inverse":       lambda v: "(1.0/%s)" % v,
               "std power":     lambda base, exp: "std::pow(%s, %s)" % (base, exp),
               "exp":           lambda v: "std::exp(%s)" % str(v),
               "ln":            lambda v: "std::log(%s)" % str(v),
               "cos":           lambda v: "std::cos(%s)" % str(v),
               "sin":           lambda v: "std::sin(%s)" % str(v),
               "absolute value":lambda v: "std::abs(%s)" % str(v),
               "sqrt":          lambda v: "std::sqrt(%s)" % str(v),
               "addition":      lambda v: _add(v),
               "multiply":      lambda v: _multiply(v),
               "power":         lambda base, exp: _power(base, exp),
               "inner product": lambda v, w: _inner_product(v, w),
               "assign":        lambda v, w: "%s = %s;" % (v, str(w)),
               "component":     lambda v, k: _component(v, k)})

# Formatting used in tabulate_tensor
format.update({"element tensor":  lambda i: "A[%d]" % i,
               "geometry tensor":
                                  lambda j, a: "G%d_%s" % (j, "_".join(["%d" % i for i in a])),
               "coefficient":     lambda j, k: format["component"]("w", [j, k]),
               "transform":       lambda t, j, k, r: _transform(t, j, k, r)})

# Geometry related variable names
format.update({"entity index": "c.entity_indices",
               "num entities": "m.num_entities",
               "cell":   lambda s: "ufc::%s" % s,
               "J":      lambda i, j: "J_%d%d" % (i, j),
               "inv(J)": lambda i, j: "K_%d%d" % (i, j),
               "det(J)": lambda r="": "detJ%s" % r})

# Code snippets
from codesnippets import *
format.update({"cell coordinates": cell_coordinates,
               "jacobian": lambda n, r="": jacobian[n] % {"restriction": r},
               "inverse jacobian": lambda n, r="": inverse_jacobian[n] % {"restriction": r},
               "jacobian and inverse": lambda n, r="": format["jacobian"](n, r) +\
                                       "\n" + format["inverse jacobian"](n, r),
               "facet determinant": lambda n, r="": facet_determinant[n] % {"restriction": r},
               "fiat coordinate map": lambda n: fiat_coordinate_map[n],
               "generate normal": lambda d, i: _generate_normal(d, i),
               "scale factor snippet": scale_factor,
               "map onto physical": map_onto_physical,
               "combinations": combinations_snippet,
               "transform snippet": transform_snippet,
               "evaluate function": evaluate_f,
               "ufc header": header_ufc,
               "dolfin header": header_dolfin,
               "footer": footer
               })

# TODO: Stuff from format_old used by KBO, should be moved around and possibly renamed.
format.update({# Loop indices
               "integration points": "ip",
               "free indices":  ["r","s","t","u"],
               "first free index": "j",
               "second free index": "k",
               "tmp value": lambda i: "tmp%d" % i,
               "tmp ref value": lambda i: "tmp_ref%d" % i,
               # Snippet variable names
               "x coordinate": "X",
               "y coordinate": "Y",
               "z coordinate": "Z",
               "scale factor": "det",
               "normal component": lambda r, j: "n%s%s" % (choose_map[r], j),
               # Random variable names
               "local dof": "dof",
               "basisvalues": "basisvalues",
               "geometry constant": "G",
               "coefficients": lambda i: "coefficients%d" %(i),
               "num derivatives": "num_derivatives",
               "derivative combinations": "combinations",
               "transform matrix": "transform",
               "transform Jinv": "Jinv",
               "dmats":   lambda i: "dmats%s" %(i),
               "dmats old": "dmats_old",
               "reference derivatives": "derivatives",
               "function value": "F",
               "nonzero columns": lambda i: "nzc%d" % i,
               "element tensor quad": "A",
               "weight": lambda i: "W%d" % (i),
               # UFC argument names
               "argument values": "values",
               "argument coordinates": "coordinates",
               "argument basis num": "i",
               "argument derivative order": "n",
               # Convenience formats
               "dof map if": lambda i,j: "%d <= %s && %s <= %d" %(i,\
                 format["argument basis num"], format["argument basis num"], j),
               # Misc
               "pointer": "*",
               "new": "new ",
#               "delete": "delete ",
               "delete pointer": lambda v, w: "delete [] %s%s;" % (v, w),
               "separator": ", ",
               "block separator": ",\n",
               "new line": "\\\n"
})

# Class names
format.update({"classname finite_element": \
                   lambda prefix, i: "%s_finite_element_%d" % (prefix.lower(), i),
               "classname dof_map": \
                   lambda prefix, i: "%s_dof_map_%d" % (prefix.lower(), i),
               "classname cell_integral": \
                   lambda prefix, form_id, sub_domain: "%s_cell_integral_%d_%d" % (prefix.lower(), form_id, sub_domain),
               "classname exterior_facet_integral": \
                   lambda prefix, form_id, sub_domain: "%s_exterior_facet_integral_%d_%d" % (prefix.lower(), form_id, sub_domain),
               "classname interior_facet_integral": \
                   lambda prefix, form_id, sub_domain: "%s_interior_facet_integral_%d_%d" % (prefix.lower(), form_id, sub_domain),
               "classname form": \
                   lambda prefix, i: "%s_form_%d" % (prefix.lower(), i)})

# Misc
format.update({"bool":    lambda v: {True: "true", False: "false"}[v],
#               "float":   lambda v: "%f" % v,
               "str":     lambda v: "%s" % str(v),
#               "epsilon": FFC_OPTIONS["epsilon"]
})

def _declaration(type, name, value=None):
    if value is None:
        return "%s %s;\n" % (type, name);
    return "%s %s = %s;\n" % (type, name, str(value));

def _component(var, k):
    if not isinstance(k, (list, tuple)):
        k = [k]
    return "%s" % var + "".join("[%s]" % str(i) for i in k)

# Utility functions for arithmetic operations
def _multiply(factors):
    """
    Generate string multiplying a list of numbers or strings.  If a
    factor is zero, the whole product is zero. Any factors equal to
    one are ignored.
    """

    # FIXME: This could probably be way more robust and elegant.

    cpp_str = format["str"]
    non_zero_factors = []
    for f in factors:

        # Round-off if f is smaller than epsilon
        if isinstance(f, (int, float)):
            if abs(f) < format["epsilon"]:
                return cpp_str(0)

        # Convert to string
        f = cpp_str(f)

        # Return zero if any factor is zero
        if f == "0":
            return cpp_str(0)

        # Ignore 1 factors
        if f == "1" or f == "1.0":
            continue

        # If sum-like, parentheseze factor
        if "+" in f or "-" in f:
            f = "(%s)" % f

        non_zero_factors += [f]

    if len(non_zero_factors) == 0:
        return cpp_str(1.0)

    return "*".join(non_zero_factors)

def _add(terms):
    "Generate string summing a list of strings."

    # FIXME: Subtract absolute value of negative numbers
    result = " + ".join([str(t) for t in terms if (str(t) != "0")])
    if result == "":
        return format["str"](0)
    return result

def _power(base, exponent):
    "Generate code for base^exponent."
    if exponent >= 0:
        return _multiply(exponent*(base,))
    else:
        return "1.0 / (%s)" % _power(base, -exponent)

def _inner_product(v, w):
    "Generate string for v[0]*w[0] + ... + v[n]*w[n]."

    # Check that v and w have same length
    assert(len(v) == len(w)), \
                  "Sizes differ (%d, %d) in inner-product!" % (len(v), len(w))

    return format["addition"]([format["multiply"]([v[i], w[i]])
                          for i in range(len(v))])

def _transform(type, j, k, r):
    map_name = {"J": "J", "JINV": "K"}[type] + {None: "", "+": "0", "-": "1"}[r]
    return (map_name + "_%d%d") % (j, k)

# FIXME: Input to _generate_switch should be a list of tuples (i, case)

def _generate_switch(variable, cases, default=None, numbers=None):
    "Generate switch statement from given variable and cases"

    # Special case: no cases and no default
    if len(cases) == 0 and default is None:
        return format["do nothing"]
    elif len(cases) == 0:
        return default

    # Special case: one case and no default
    if len(cases) == 1 and default is None:
        return cases[0]

    # Create numbers for switch
    if numbers is None:
        numbers = range(len(cases))

    # Create switch
    code = "switch (%s)\n{\n" % variable
    for (i, case) in enumerate(cases):
        code += "case %d:\n  {\n  %s\n    break;\n  }\n" % (numbers[i], indent(case, 2))
    code += "}\n"

    # Default value
    if default:
        code += "\n" + default

    return code

def inner_product(a, b, format):
    """Generate code for inner product of a and b, where a is a list
    of floating point numbers and b is a list of symbols."""

    # Check input
    if not len(a) == len(b):
        error("Dimensions don't match for inner product.")

    # Prefetch formats to speed up code generation
    format_add            = format["addition"]
    format_subtract       = format["subtract"]
    format_multiply       = format["multiply"]
    format_floating_point = format["floating point"]
    format_epsilon        = format["epsilon"]

    # Add all entries
    value = None
    for i in range(len(a)):

        # Skip terms where a is almost zero
        if abs(a[i]) <= format_epsilon:
            continue

        # Fancy handling of +, -, +1, -1
        if value:
            if abs(a[i] - 1.0) < format_epsilon:
                value = format_add([value, b[i]])
            elif abs(a[i] + 1.0) < format_epsilon:
                value = format_subtract([value, b[i]])
            elif a[i] > 0.0:
                value = format_add([value, format_multiply([format_floating_point(a[i]), b[i]])])
            else:
                value = format_subtract([value, format_multiply([format_floating_point(-a[i]), b[i]])])
        else:
            if abs(a[i] - 1.0) < format_epsilon or abs(a[i] + 1.0) < format_epsilon:
                value = b[i]
            else:
                value = format_multiply([format_floating_point(a[i]), b[i]])

    return value or format_floating_point(0.0)


def tabulate_matrix(matrix, format):
    "Function that tabulates the values of a matrix, into a two dimensional array."

    # Check input
    if not len(numpy.shape(matrix)) == 2:
        error("This is not a matrix.")

    # Prefetch formats to speed up code generation
    format_block          = format["block"]
    format_separator      = format["separator"]
    format_floating_point = format["floating point"]
    format_epsilon        = format["epsilon"]

    # Get size of matrix
    num_rows = numpy.shape(matrix)[0]
    num_cols = numpy.shape(matrix)[1]

    # Set matrix entries equal to zero if their absolute values is smaller than format_epsilon
    for i in range(num_rows):
        for j in range(num_cols):
            if abs(matrix[i][j]) < format_epsilon:
                matrix[i][j] = 0.0

    # Generate array of values
    value = format["new line"] + format["block begin"]
    rows = []

    for i in range(num_rows):
        rows += [format_block(format_separator.join([format_floating_point(matrix[i,j])\
                 for j in range(num_cols)]))]

    value += format["block separator"].join(rows)
    value += format["block end"]

    return value

def tabulate_vector(vector, format):
    "Function that tabulates the values of a vector, into a one dimensional array."

    # Check input
    if not len(numpy.shape(vector)) == 1:
        error("This is not a vector.")

    # Prefetch formats to speed up code generation
    format_block          = format["block"]
    format_separator      = format["separator"]
    format_floating_point = format["floating point"]
    format_epsilon        = format["epsilon"]

    # Get size of matrix
    num_cols = numpy.shape(vector)[0]

    # Set vector entries equal to zero if their absolute values is smaller than format_epsilon
    for i in range(num_cols):
        if abs(vector[i]) < format_epsilon:
            vector[i] = 0.0

    value = format_block(format_separator.join([format_floating_point(val) for val in vector]))

    return value


# ---- Indentation control ----
class IndentControl:
    "Class to control the indentation of code"

    def __init__(self):
        "Constructor"
        self.size = 0
        self.increment = 2

    def increase(self):
        "Increase indentation by increment"
        self.size += self.increment

    def decrease(self):
        "Decrease indentation by increment"
        self.size -= self.increment

    def indent(self, a):
        "Indent string input string by size"
        return indent(a, self.size)

def indent(block, num_spaces):
    "Indent each row of the given string block with n spaces."
    indentation = " " * num_spaces
    return indentation + ("\n" + indentation).join(block.split("\n"))



# FIXME: Major cleanup needed, remove as much as possible
#from codesnippets import *

# FIXME: KBO: temporary hack to get dictionary working.
from constants import FFC_OPTIONS
import platform
options=FFC_OPTIONS.copy()

# Old dictionary, move the stuff we need to the new dictionary above
choose_map = {None: "", "+": "0", "-": 1}
transform_options = {"JINV": lambda m, j, k: "Jinv%s_%d%d" % (m, j, k),
                     "J": lambda m, j, k: "J%s_%d%d" % (m, k, j)}
format_old = {
    # Operators
    #
#    "multiply": lambda v: _multiply(v),
    "multiply": lambda v: "*".join(v),
    "times equal": lambda i, j: "%s *= %s;" %(i, j),
    "add equal": lambda i, j: "%s += %s;" % (i, j),
    "inverse": lambda v: "(1.0/%s)" % v,
    "absolute value": lambda v: "std::abs(%s)" % v,
    "sqrt": lambda v: "std::sqrt(%s)" % v,
    "add": lambda v: " + ".join(v),
    "subtract": lambda v: " - ".join(v),
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
    "grouping": lambda v: "(%s)" % v,
    "block": lambda v: "{%s}" % v,
    "block begin": "{",
    "block end": "}",
    "separator": ", ",
    "block separator": ",\n",
    "generate jacobian": lambda d, i: _generate_jacobian(d, i),
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
#    "const uint declaration": "const unsigned int ",
    "static const uint declaration": "static const unsigned int ",
    "static uint declaration": "static unsigned int ",
    "table declaration": "static const double ",
    # variable names
#    "element tensor quad": "A",
#    "integration points": "ip",
#    "first free index": "j",
#    "second free index": "k",
#    "free secondary indices":["r","s","t","u"],
    "derivatives": lambda i,j,k,l: "dNdx%d_%d[%s][%s]" % (i,j,k,l),
    "element coordinates": lambda i,j: "x[%s][%s]" % (i,j),
    "weight": lambda i: "W%d" % (i),
    "weights": lambda i,j: self.format["weight"](i) + "[%s]" % (j),
    "psis": "P",
#    "function value": "F",
#    "argument coordinates": "coordinates",
#    "argument values": "values",
#    "argument basis num": "i",
#    "argument derivative order": "n",
#    "local dof": "dof",
#    "x coordinate": "x",
#    "y coordinate": "y",
#    "z coordinate": "z",
#    "fiat x coordinate": "fiat_x",
#    "fiat y coordinate": "fiat_y",
#    "fiat z coordinate": "fiat_z",
#    "scalings": lambda i,j: "scalings_%s_%d" %(i,j),
#    "coefficients table": lambda i: "coefficients%d" %(i),
#    "basisvalues table": "basisvalues",
#    "basisvalues": lambda i: "basisvalues"  + format_old["array access"](i),
#    "dmats table": lambda i: "dmats%s" %(i),
#    "dmats": "dmats",
#    "dmats old": "dmats_old",
#    "coefficient scalar": lambda i: "coeff%d" %(i),
#    "new coefficient scalar": lambda i: "new_coeff%d" %(i),
#    "psitilde_a": "psitilde_a",
#    "psitilde_bs": lambda i: "psitilde_bs_%d" %(i),
#    "psitilde_cs": lambda i,j: "psitilde_cs_%d%d" %(i,j),
#    "basisvalue": lambda i: "basisvalue%d" %(i),
#    "evaluate_basis aux index": lambda i: "idx%d" %(i),
#    "evaluate_basis aux factor": lambda i: "f%d" %(i),
#    "evaluate_basis aux value": lambda i: "n%d" %(i),
#    "num derivatives": "num_derivatives",
#    "reference derivatives": "derivatives",
#    "derivative combinations": "combinations",
#    "transform matrix": "transform",
#    "transform Jinv": "Jinv",
    "normal component": lambda r, j: "n%s%s" % (choose_map[r], j),
    "tmp declaration": lambda j, k: "const double " + format_old["tmp access"](j, k),
    "tmp access": lambda j, k: "tmp%d_%d" % (j, k),
    "determinant": lambda r: "detJ%s" % choose_map[r],
    "scale factor": "det",
#    "constant": lambda j: "c%d" % j,
#    "coefficient table": lambda j, k: "w[%d][%d]" % (j, k),
#    "coefficient": lambda j, k: "w[%d][%d]" % (j, k),
#    "coeff": "w",
#    "modified coefficient declaration": lambda i, j, k, l: "const double c%d_%d_%d_%d" % (i, j, k, l),
#    "modified coefficient access": lambda i, j, k, l: "c%d_%d_%d_%d" % (i, j, k, l),
#    "transform": lambda type, j, k, r: "%s" % (transform_options[type](choose_map[r], j, k)),
#    "transform_ufl": lambda type, j, k, r: "%s" % (transform_options_ufl[type](choose_map[r], j, k)),
#    "reference tensor" : lambda j, i, a: None,
#    "geometry tensor": "G",
#    "sign tensor": lambda type, i, k: "S%s%s_%d" % (type, i, k),
#    "sign tensor declaration": lambda s: "const int " + s,
#    "signs": "S",
#    "vertex values": lambda i: "vertex_values[%d]" % i,
#    "dof values": lambda i: "dof_values[%d]" % i,
#    "dofs": lambda i: "dofs[%d]" % i,
#    "entity index": lambda d, i: "c.entity_indices[%d][%d]" % (d, i),
#    "num entities": lambda dim : "m.num_entities[%d]" % dim,
#    "offset declaration": "unsigned int offset",
#    "offset access": "offset",
#    "nonzero columns": lambda i: "nzc%d" % i,
    # access
#    "array access": lambda i: "[%s]" %(i),
#    "matrix access": lambda i,j: "[%s][%s]" %(i,j),
    "secondary index": lambda i: "_%s" %(i),
    # program flow
    "dof map if": lambda i,j: "if (%d <= %s && %s <= %d)" %(i,\
                                                                format_old["argument basis num"], format_old["argument basis num"], j),
    "loop": lambda i,j,k: "for (unsigned int %s = %s; %s < %s; %s++)"% (i, j, i, k, i),
    "if": "if",
    # snippets
    "coordinate map": lambda s: eval("map_coordinates_%s" % s),
    "coordinate map FIAT": lambda s: eval("map_coordinates_FIAT_%s" % s),
    "facet sign": lambda e: "sign_facet%d" % e,
#    "snippet facet signs": lambda d: eval("facet_sign_snippet_%dD" % d),
#    "snippet dof map": evaluate_basis_dof_map,
#    "snippet eta_interval": eta_interval_snippet,
#    "snippet eta_triangle": eta_triangle_snippet,
#    "snippet eta_tetrahedron": eta_tetrahedron_snippet,
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
    #"snippet calculate dof": calculate_dof,
    "get cell vertices" : "const double * const * x = c.coordinates;",
    "generate body": lambda d: _generate_body(d),
    # misc
    "comment": lambda v: "// %s" % v,
#    "pointer": "*",
#    "new": "new ",
#    "delete": "delete ",
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
#precision = int(options["precision"])
#f1 = "%%.%dg" % precision
#f2 = "%%.%de" % precision

#def floating_point(v):
#    "Format floating point number."
#    if abs(v) < 100.0:
#        return f1 % v
#    else:
#        return f2 % v

#def floating_point_windows(v):
#    "Format floating point number for Windows (remove extra leading zero in exponents)."
#    return floating_point(v).replace("e-0", "e-").replace("e+0", "e+")

#if platform.system() == "Windows":
#    format["floating point"] = floating_point_windows
#    format_old["floating point"] = floating_point_windows
#else:
#    format_old["floating point"] = floating_point
#    format["floating point"] = floating_point

#format_old["epsilon"] = 10.0*eval("1e-%s" % precision)
#format["epsilon"] = 10.0*eval("1e-%s" % precision)

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

# Declarations to examine
types = [["double"],
         ["const", "double"],
         ["const", "double", "*", "const", "*"],
         ["int"],
         ["const", "int"],
         ["unsigned", "int"],
         ["bool"],
         ["const", "bool"]]

# Special characters and delimiters
special_characters = ["+", "-", "*", "/", "=", ".", " ", ";", "(", ")", "\\", "{", "}", "[","]"]

def remove_unused(code, used_set=set()):
    """
    Remove unused variables from a given C++ code. This is useful when
    generating code that will be compiled with gcc and options -Wall
    -Werror, in which case gcc returns an error when seeing a variable
    declaration for a variable that is never used.

    Optionally, a set may be specified to indicate a set of variables
    names that are known to be used a priori.
    """

    # Dictionary of (declaration_line, used_lines) for variables
    variables = {}

    # List of variable names (so we can search them in order)
    variable_names = []

    lines = code.split("\n")
    for (line_number, line) in enumerate(lines):
        # Split words
        words = [word for word in line.split(" ") if not word == ""]
        # Remember line where variable is declared
        for type in [type for type in types if " ".join(type) in " ".join(words)]: # Fewer matches than line below.
        # for type in [type for type in types if len(words) > len(type)]:
            variable_type = words[0:len(type)]
            variable_name = words[len(type)]

            # Skip special characters
            if variable_name in special_characters:
                continue

            # Test if any of the special characters are present in the variable name
            # If this is the case, then remove these by assuming that the 'real' name
            # is the first entry in the return list. This is implemented to prevent
            # removal of e.g. 'double array[6]' if it is later used in a loop as 'array[i]'
            if variable_type == type:

                # Create correct variable name (e.g. y instead of
                # y[2]) for variables with separators
                seps_present = [sep for sep in special_characters if sep in variable_name]
                if seps_present:
                    variable_name = [variable_name.split(sep)[0] for sep in seps_present]
                    variable_name.sort()
                    variable_name = variable_name[0]

                variables[variable_name] = (line_number, [])
                if not variable_name in variable_names:
                    variable_names += [variable_name]

        # Mark line for used variables
        for variable_name in variables:
            (declaration_line, used_lines) = variables[variable_name]
            if _variable_in_line(variable_name, line) and line_number > declaration_line:
                variables[variable_name] = (declaration_line, used_lines + [line_number])

    # Reverse the order of the variable names to catch variables used
    # only by variables that are removed
    variable_names.reverse()

    # Remove declarations that are not used
    removed_lines = []
    for variable_name in variable_names:
        (declaration_line, used_lines) = variables[variable_name]
        for line in removed_lines:
            if line in used_lines:
                used_lines.remove(line)
        if not used_lines and not variable_name in used_set:
            debug("Removing unused variable: %s" % variable_name)
            lines[declaration_line] = None # KBO: Need to completely remove line for evaluate_basis* to work
            # lines[declaration_line] = "// " + lines[declaration_line]
            removed_lines += [declaration_line]
    return "\n".join([line for line in lines if not line is None])

def count_ops(code):
    "Count the number of operations in code (multiply-add pairs)."
    num_add = code.count("+") + code.count("-")
    num_multiply = code.count("*") + code.count("/")
    return (num_add + num_multiply) / 2

def _variable_in_line(variable_name, line):
    "Check if variable name is used in line"
    if not variable_name in line:
        return False
    for character in special_characters:
        line = line.replace(character, "\\" + character)
    delimiter = "[" + ",".join(["\\" + c for c in special_characters]) + "]"
    return not re.search(delimiter + variable_name + delimiter, line) == None

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

def _generate_normal(geometric_dimension, domain_type, reference_normal=False):
    "Generate code for computing normal"

    # Choose snippets
    direction = normal_direction[geometric_dimension]
    normal = facet_normal[geometric_dimension]

    # Choose restrictions
    if domain_type == "exterior_facet":
        code = direction % {"restriction": "", "facet" : "facet"}
        code += normal % {"direction" : "", "restriction": ""}
    elif domain_type == "interior_facet":
        code = direction % {"restriction": choose_map["+"], "facet": "facet0"}
        code += normal % {"direction" : "", "restriction": choose_map["+"]}
        code += normal % {"direction" : "!", "restriction": choose_map["-"]}
    else:
        error("Unsupported domain_type: %s" % str(domain_type))
    return code

def set_float_formatting(precision):
    "Set floating point formatting based on precision."

    # Options for float formatting
    f1 = "%%.%df" % precision
    f2 = "%%.%de" % precision

    # Regular float formatting
    def floating_point_regular(v):
        if abs(v) < 100.0:
            return f1 % v
        else:
            return f2 % v

    # Special float formatting on Windows (remove extra leading zero)
    def floating_point_windows(v):
        return floating_point(v).replace("e-0", "e-").replace("e+0", "e+")

    # Set float formatting
    if platform.system() == "Windows":
        format["float"] = floating_point_windows
    else:
        format["float"] = floating_point_regular

    # FIXME: KBO: Remove once we agree on the format of 'f1'
    format["floating point"] = format["float"]

    # Set machine precision
    format["epsilon"] = 10.0*eval("1e-%s" % precision)

def set_exception_handling(convert_exceptions_to_warnings):
    "Set handling of exceptions."
    if convert_exceptions_to_warnings:
        format["exception"] = format["warning"]
