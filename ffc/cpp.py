"This module defines rules and algorithms for generating C++ code."

__author__ = "Anders Logg (logg@simula.no) and friends"
__date__ = "2009-12-16"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard 2010
# Modified by Marie E. Rognes 2010
# Last changed: 2010-03-11

# Python modules
import re, numpy, platform

# FFC modules
from ffc.log import debug, error

# Formatting rules
# FIXME: KBO: format is a builtin_function, i.e., we should use a different name.
format = {}
choose_map = {None: "", "+": "0", "-": "1"}

# Program flow
format.update({
    "return":         lambda v: "return %s;" % str(v),
    "grouping":       lambda v: "(%s)" % v,
    "block":          lambda v: "{%s}" % v,
    "block begin":    "{",
    "block end":      "}",
    "list":           lambda v: format["block"](format["list separator"].join([str(l) for l in v])),
    "switch":         lambda v, cases, default=None, numbers=None: _generate_switch(v, cases, default, numbers),
    "exception":      lambda v: "throw std::runtime_error(\"%s\");" % v,
    "warning":        lambda v: 'std::cerr << "*** FFC warning: " << "%s" << std::endl;' % v,
    "comment":        lambda v: "// %s" % v,
    "if":             lambda c, v: "if (%s)\n{\n%s\n}\n" % (c, v),
    "loop":           lambda i, j, k: "for (unsigned int %s = %s; %s < %s; %s++)"% (i, j, i, k, i),
    "generate loop":  lambda v, w, _indent=0: _generate_loop(v, w, _indent),
    "is equal":       " == ",
    "not equal":      " == ",
    "less than":      " < ",
    "greater than":   " > ",
    "less equal":     " <= ",
    "greater equal":  " >= ",
    "do nothing":     "// Do nothing"
})

# Declarations
format.update({
    "declaration":                    lambda t, n, v=None: _declaration(t, n, v),
    "float declaration":              "double",
    "int declaration":                "int",
    "uint declaration":               "unsigned int",
    "static const uint declaration":  "static const unsigned int",
    "static const float declaration": "static const double",
    "const float declaration":        lambda v, w: "const double %s = %s;" % (v, w),
    "const uint declaration":         lambda v, w: "const unsigned int %s = %s;" % (v, w),
    "dynamic array":                  lambda t, n, s: "%s *%s = new %s[%s];" % (t, n, t, s),
    "delete dynamic array":           lambda n, s=None: _delete_array(n, s),
    "create foo":                     lambda v: "new %s" % v
})

# Mathematical operators
format.update({
    "add":            lambda v: " + ".join(v),
    "iadd":           lambda v, w: "%s += %s;" % (str(v), str(w)),
    "sub":            lambda v: " - ".join(v),
    "neg":            lambda v: "-%s" % v,
    "mul":            lambda v: "*".join(v),
    "imul":           lambda v, w: "%s *= %s;" % (str(v), str(w)),
    "div":            lambda v, w: "%s/%s" % (str(v), str(w)),
    "inverse":        lambda v: "(1.0/%s)" % v,
    "std power":      lambda base, exp: "std::pow(%s, %s)" % (base, exp),
    "exp":            lambda v: "std::exp(%s)" % str(v),
    "ln":             lambda v: "std::log(%s)" % str(v),
    "cos":            lambda v: "std::cos(%s)" % str(v),
    "sin":            lambda v: "std::sin(%s)" % str(v),
    "tan":            lambda v: "std::tan(%s)" % str(v),
    "acos":           lambda v: "std::acos(%s)" % str(v),
    "asin":           lambda v: "std::asin(%s)" % str(v),
    "atan":           lambda v: "std::atan(%s)" % str(v),
    "absolute value": lambda v: "std::abs(%s)" % str(v),
    "sqrt":           lambda v: "std::sqrt(%s)" % str(v),
    "addition":       lambda v: _add(v),
    "multiply":       lambda v: _multiply(v),
    "power":          lambda base, exp: _power(base, exp),
    "inner product":  lambda v, w: _inner_product(v, w),
    "assign":         lambda v, w: "%s = %s;" % (v, str(w)),
    "component":      lambda v, k: _component(v, k)
})

# Formatting used in tabulate_tensor
format.update({
    "geometry tensor": lambda j, a: "G%d_%s" % (j, "_".join(["%d" % i for i in a]))
})

# Geometry related variable names (from code snippets).
format.update({
    "entity index":     "c.entity_indices",
    "num entities":     "m.num_entities",
    "cell":             lambda s: "ufc::%s" % s,
    "J":                lambda i, j: "J_%d%d" % (i, j),
    "inv(J)":           lambda i, j: "K_%d%d" % (i, j),
    "det(J)":           lambda r=None: "detJ%s" % choose_map[r],
    "cell volume":      lambda r=None: "volume%s" % choose_map[r],
    "circumradius":     lambda r=None: "circumradius%s" % choose_map[r],
    "scale factor":     "det",
    "transform":        lambda t, j, k, r: _transform(t, j, k, r),
    "normal component": lambda r, j: "n%s%s" % (choose_map[r], j),
    "x coordinate":     "X",
    "y coordinate":     "Y",
    "z coordinate":     "Z",
    "ip coordinates":   lambda i, j: "X%d[%d]" % (i, j),
    "affine map table": lambda i, j: "FEA%d_f%d" % (i, j),
    "coordinates":      "x"
})

# UFC function arguments and class members (names)
format.update({
    "element tensor":             lambda i: "A[%s]" % i,
    "element tensor term":        lambda i, j: "A%d[%s]" % (j, i),
    "coefficient":                lambda j, k: format["component"]("w", [j, k]),
    "argument basis num":         "i",
    "argument derivative order":  "n",
    "argument values":            "values",
    "argument coordinates":       "coordinates",
    "facet":                      lambda r: "facet%s" % choose_map[r],
    "argument axis":              "i",
    "argument dimension":         "d",
    "argument entity":            "i",
    "member global dimension":    "_global_dimension",
    "argument dofs":              "dofs",
    "argument dof num":           "i",
    "argument dof values":        "dof_values",
    "argument vertex values":     "vertex_values",
    "argument sub":               "i" # sub domain, sub element
})

# Formatting used in evaluatedof.
format.update({
    "dof vals":                 "vals",
    "dof result":               "result",
    "dof X":                    lambda i: "X_%d" % i,
    "dof D":                    lambda i: "D_%d" % i,
    "dof W":                    lambda i: "W_%d" % i,
    "dof copy":                 lambda i: "copy_%d" % i,
    "dof physical coordinates": "y"
})


# Formatting used in evaluate_basis, evaluate_basis_derivatives and quadrature
# code generators.
format.update({
    # evaluate_basis and evaluate_basis_derivatives
    "tmp value":                lambda i: "tmp%d" % i,
    "tmp ref value":            lambda i: "tmp_ref%d" % i,
    "local dof":                "dof",
    "basisvalues":              "basisvalues",
    "coefficients":             lambda i: "coefficients%d" %(i),
    "num derivatives":          "num_derivatives",
    "derivative combinations":  "combinations",
    "transform matrix":         "transform",
    "transform Jinv":           "Jinv",
    "dmats":                    lambda i: "dmats%s" %(i),
    "dmats old":                "dmats_old",
    "reference derivatives":    "derivatives",
    "dof values":               "dof_values",
    "dof map if":               lambda i,j: "%d <= %s && %s <= %d"\
                                % (i, format["argument basis num"], format["argument basis num"], j),
    "dereference pointer":      lambda n: "*%s" % n,
    "reference variable":       lambda n: "&%s" % n,
    "call basis":               lambda i, s: "evaluate_basis(%s, %s, coordinates, c);" % (i, s),
    "call basis_derivatives":   lambda i, s: "evaluate_basis_derivatives(%s, n, %s, coordinates, c);" % (i, s),

    # quadrature code generators
    "integration points": "ip",
    "first free index":   "j",
    "second free index":  "k",
    "geometry constant":  lambda i: "G[%d]" % i,
    "ip constant":        lambda i: "I[%d]" % i,
    "basis constant":     lambda i: "B[%d]" % i,
    "conditional":        lambda i: "C[%d]" % i,
    "evaluate conditional":lambda i,j,k: "(%s) ? %s : %s" % (i,j,k),
#    "geometry constant":  lambda i: "G%d" % i,
#    "ip constant":        lambda i: "I%d" % i,
#    "basis constant":     lambda i: "B%d" % i,
    "function value":     lambda i: "F%d" % i,
    "nonzero columns":    lambda i: "nzc%d" % i,
    "weight":             lambda i: "W%d" % (i),
    "psi name":           lambda c, f, co, d: _generate_psi_name(c,f,co,d),
    # both
    "free indices":       ["r","s","t","u"],
    "matrix index":       lambda i, j, range_j: _matrix_index(i, str(j), str(range_j))
})

# Misc
format.update({
    "bool":             lambda v: {True: "true", False: "false"}[v],
    "str":              lambda v: "%s" % v,
    "int":              lambda v: "%d" % v,
    "list separator":   ", ",
    "block separator":  ",\n",
    "new line":         "\\\n",
    "tabulate tensor":  lambda m: _tabulate_tensor(m),
})

# Code snippets
from codesnippets import *
format.update({
    "cell coordinates":     cell_coordinates,
    "jacobian":             lambda n, r="": jacobian[n] % {"restriction": r},
    "inverse jacobian":     lambda n, r="": inverse_jacobian[n] % {"restriction": r},
    "jacobian and inverse": lambda n, r=None: format["jacobian"](n, choose_map[r]) +\
                            "\n" + format["inverse jacobian"](n, choose_map[r]),
    "facet determinant":    lambda n, r=None: facet_determinant[n] % {"restriction": choose_map[r]},
    "fiat coordinate map":  lambda n: fiat_coordinate_map[n],
    "generate normal":      lambda d, i: _generate_normal(d, i),
    "generate cell volume": lambda d, i: _generate_cell_volume(d, i),
    "generate circumradius": lambda d, i: _generate_circumradius(d, i),
    "generate ip coordinates":  lambda g, num_ip, name, ip, r=None: (ip_coordinates[g][0], ip_coordinates[g][1] % \
                                {"restriction": choose_map[r], "ip": ip, "name": name, "num_ip": num_ip}),
    "scale factor snippet": scale_factor,
    "map onto physical":    map_onto_physical,
    "combinations":         combinations_snippet,
    "transform snippet":    transform_snippet,
    "evaluate function":    evaluate_f,
    "ufc comment":          comment_ufc,
    "dolfin comment":       comment_dolfin,
    "header_h":             header_h,
    "header_c":             header_c,
    "footer":               footer
})

# Class names
format.update({
    "classname finite_element": lambda prefix, i:\
               "%s_finite_element_%d" % (prefix.lower(), i),

    "classname dof_map":  lambda prefix, i: "%s_dof_map_%d" % (prefix.lower(), i),

    "classname cell_integral":  lambda prefix, form_id, sub_domain:\
               "%s_cell_integral_%d_%d" % (prefix.lower(), form_id, sub_domain),

    "classname exterior_facet_integral":  lambda prefix, form_id, sub_domain:\
              "%s_exterior_facet_integral_%d_%d" % (prefix.lower(), form_id, sub_domain),

    "classname interior_facet_integral":  lambda prefix, form_id, sub_domain:\
              "%s_interior_facet_integral_%d_%d" % (prefix.lower(), form_id, sub_domain),

    "classname form": lambda prefix, i: "%s_form_%d" % (prefix.lower(), i)
})

# Helper functions for formatting.
def _declaration(type, name, value=None):
    if value is None:
        return "%s %s;" % (type, name);
    return "%s %s = %s;" % (type, name, str(value));

def _component(var, k):
    if not isinstance(k, (list, tuple)):
        k = [k]
    return "%s" % var + "".join("[%s]" % str(i) for i in k)

def _delete_array(name, size=None):
    if size is None:
        return "delete [] %s;" % name
    f_r = format["free indices"][0]
    code = format["generate loop"](["delete [] %s;" % format["component"](name, f_r)], [(f_r, 0, size)])
    code.append("delete [] %s;" % name)
    return "\n".join(code)

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
            if abs(f - 1.0) < format["epsilon"]:
                continue

        # Convert to string
        f = cpp_str(f)

        # Return zero if any factor is zero
        if f == "0":
            return cpp_str(0)

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
    assert(len(v) == len(w)), "Sizes differ in inner-product!"

    # Special case, zero terms
    if len(v) == 0: return format["float"](0)

    # Straightforward handling when we only have strings
    if isinstance(v[0], str):
        return _add([_multiply([v[i], w[i]]) for i in range(len(v))])

    # Fancy handling of negative numbers etc
    result = None
    eps = format["epsilon"]
    add = format["add"]
    sub = format["sub"]
    neg = format["neg"]
    mul = format["mul"]
    fl  = format["float"]
    for (c, x) in zip(v, w):
        if result:
            if abs(c - 1.0) < eps:
                result = add([result, x])
            elif abs(c + 1.0) < eps:
                result = sub([result, x])
            elif c > eps:
                result = add([result, mul([fl(c), x])])
            elif c < -eps:
                result = sub([result, mul([fl(-c), x])])
        else:
            if abs(c - 1.0) < eps:
                result = x
            elif abs(c + 1.0) < eps:
                result = neg(x)
            elif c > eps:
                result = mul([fl(c), x])
            elif c < -eps:
                result = neg(mul([fl(-c), x]))

    return result

def _transform(type, j, k, r):
    map_name = {"J": "J", "JINV": "K"}[type] + choose_map[r]
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

def _tabulate_tensor(vals):
    "Tabulate a multidimensional tensor. (Replace tabulate_matrix and tabulate_vector)."

    # Prefetch formats to speed up code generation
    f_block     = format["block"]
    f_list_sep  = format["list separator"]
    f_block_sep = format["block separator"]
    # FIXME: KBO: Change this to "float" once issue in set_float_formatting is fixed.
    f_float     = format["floating point"]
    f_epsilon   = format["epsilon"]

    # Create numpy array and get shape.
    tensor = numpy.array(vals)
    shape = numpy.shape(tensor)
    if len(shape) == 1:
        # Create zeros if value is smaller than tolerance.
        values = []
        for v in tensor:
            if abs(v) < f_epsilon:
                values.append(f_float(0.0))
            else:
                values.append(f_float(v))
        # Format values.
        return f_block(f_list_sep.join(values))
    elif len(shape) > 1:
        return f_block(f_block_sep.join([_tabulate_tensor(tensor[i]) for i in range(shape[0])]))
    else:
        error("Not an N-dimensional array:\n%s" % tensor)

def _generate_loop(lines, loop_vars, _indent):
    "This function generates a loop over a vector or matrix."

    # Prefetch formats to speed up code generation.
    f_loop     = format["loop"]
    f_begin    = format["block begin"]
    f_end      = format["block end"]
    f_comment  = format["comment"]

    if not loop_vars:
        return lines

    code = []
    for ls in loop_vars:
        # Get index and lower and upper bounds.
        index, lower, upper = ls
        # Loop index.
        code.append(indent(f_loop(index, lower, upper), _indent))
        code.append(indent(f_begin, _indent))

        # Increase indentation.
        _indent += 2

        # If this is the last loop, write values.
        if index == loop_vars[-1][0]:
            for l in lines:
                code.append(indent(l, _indent))

    # Decrease indentation and write end blocks.
    indices = [var[0] for var in loop_vars]
    indices.reverse()
    for index in indices:
        _indent -= 2
        code.append(indent(f_end + f_comment("end loop over '%s'" % index), _indent))

    return code

def _matrix_index(i, j, range_j):
    "Map the indices in a matrix to an index in an array i.e., m[i][j] -> a[i*range(j)+j]"
    if i == 0:
        access = j
    elif i == 1:
        access = format["add"]([range_j, j])
    else:
        irj = format["mul"]([format["str"](i), range_j])
        access = format["add"]([irj, j])
    return access

def _generate_psi_name(counter, facet, component, derivatives):
    """Generate a name for the psi table of the form:
    FE#_f#_C#_D###, where '#' will be an integer value.

    FE  - is a simple counter to distinguish the various bases, it will be
          assigned in an arbitrary fashion.

    f   - denotes facets if applicable, range(element.num_facets()).

    C   - is the component number if any (this does not yet take into account
          tensor valued functions)

    D   - is the number of derivatives in each spatial direction if any. If the
          element is defined in 3D, then D012 means d^3(*)/dydz^2."""

    name = "FE%d" % counter
    if not facet is None:
        name += "_f%d" % facet
    if component != () and component != []:
        name += "_C%d" % component
    if any(derivatives):
        name += "_D" + "".join([str(d) for d in derivatives])

    return name

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

def _generate_cell_volume(geometric_dimension, domain_type):
    "Generate code for computing cell volume."

    # Choose snippets
    volume = cell_volume[geometric_dimension]

    # Choose restrictions
    if domain_type in ("cell", "exterior_facet"):
        code = volume % {"restriction": ""}
    elif domain_type == "interior_facet":
        code = volume % {"restriction": choose_map["+"]}
        code += volume % {"restriction": choose_map["-"]}
    else:
        error("Unsupported domain_type: %s" % str(domain_type))
    return code

def _generate_circumradius(geometric_dimension, domain_type):
    "Generate code for computing a cell's circumradius."

    # Choose snippets
    radius = circumradius[geometric_dimension]

    # Choose restrictions
    if domain_type in ("cell", "exterior_facet"):
        code = radius % {"restriction": ""}
    elif domain_type == "interior_facet":
        code = radius % {"restriction": choose_map["+"]}
        code += radius % {"restriction": choose_map["-"]}
    else:
        error("Unsupported domain_type: %s" % str(domain_type))
    return code

# Functions.
def indent(block, num_spaces):
    "Indent each row of the given string block with n spaces."
    indentation = " " * num_spaces
    return indentation + ("\n" + indentation).join(block.split("\n"))

def count_ops(code):
    "Count the number of operations in code (multiply-add pairs)."
    num_add = code.count(" + ") + code.count(" - ")
    num_multiply = code.count("*") + code.count("/")
    return (num_add + num_multiply) / 2

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
    generating code that will be compiled with gcc and parameters -Wall
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
        # Exclude commented lines.
        if line[:2] == "//" or line[:3] == "///":
            continue

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

def _variable_in_line(variable_name, line):
    "Check if variable name is used in line"
    if not variable_name in line:
        return False
    for character in special_characters:
        line = line.replace(character, "\\" + character)
    delimiter = "[" + ",".join(["\\" + c for c in special_characters]) + "]"
    return not re.search(delimiter + variable_name + delimiter, line) == None
