
import os
from uflacs.utils.log import info
from uflacs.utils.str_utils import format_list, format_dict

from uflacs.codeutils.format_code import format_code, Block, Indented, Namespace

def add_compile_dolfin_options(opts):
    "Args: list of .ufl file(s)."
    pass

class DolfinExpressionFormatter(object):
    """Cpp formatter class for dolfin::Expression generation."""

    def spatial_coordinate(self, o, component=(), derivatives=()):
        if len(derivatives) > 1:
            return "0"
        i = component[0] if component else 0
        if derivatives:
            d, = derivatives
            return "1" if i == d else "0"
        else:
            return "x[%d]" % i

    # FIXME: Implement rules for geometry (facet normal, cell volume, cell circumradius)

    def _facet_normal(self, o, component=(), derivatives=()):
        if derivatives:
            return "0"
        i = component[0] if component else 0
        return "n[%d]" % i

    def _cell_volume(self, o, component=(), derivatives=()):
        uflacs_assert(not component, "Expecting no component for scalar value.")
        if derivatives:
            return "0"
        return "K_vol"

    def _circumradius(self, o, component=(), derivatives=()):
        uflacs_assert(not component, "Expecting no component for scalar value.")
        if derivatives:
            return "0"
        return "K_rad"

    # FIXME: Implement rules for coefficients

    def _coefficient(self, o, component=(), derivatives=()):
        basename = "dw" if derivatives else "w"
        code = "%s[%d]" % (basename, o.count())
        if component:
            code += " " + "".join("[%d]" % c for c in component)
        if derivatives:
            uflacs_assert(len(derivatives) == 1, "Not implemented.")
            code += " " + "".join("[%d]" % d for d in derivatives)
        return code

    def _argument(self, o, component=(), derivatives=()):
        basename = "dv" if derivatives else "v"
        code = "%s[%d]" % (basename, o.count())
        if component:
            code += " " + "".join("[%d]" % c for c in component)
        if derivatives:
            uflacs_assert(len(derivatives) == 1, "Not implemented.")
            code += " " + "".join("[%d]" % d for d in derivatives)
        return code

def compile_dolfin_expression_body(expr):
    from ufl.classes import Terminal, Indexed, SpatialDerivative
    from ufl.algorithms import Graph, preprocess_expression
    from uflacs.codeutils.cpp_format import CppFormatterRules
    from uflacs.codeutils.code_formatter import CodeFormatter

    # Construct a specialized C++ expression formatter
    target_formatter = DolfinExpressionFormatter()
    cpp_formatter = CppFormatterRules(target_formatter)
    variables = {}
    code_formatter = CodeFormatter(cpp_formatter, variables)

    # Preprocess expression and build computational graph
    expr_data = preprocess_expression(expr)
    expr = expr_data.preprocessed_expr
    G = Graph(expr)
    V, E = G

    # Generate code for graph vertices
    listing = []
    k = 0
    for i, v in enumerate(V):
        if not isinstance(v, (Terminal, Indexed, SpatialDerivative)):
            vname = 's[%d]' % k
            k += 1
            vcode = code_formatter.visit(v)
            listing.append('%s = %s; // %s' % (vname, vcode, str(v)))
            variables[v] = vname

    # Generate code prelude
    # FIXME: Compute coefficients
    # FIXME: Compute geometry
    prelude = ['double s[%d];' % k]

    # Generate assignments to values[] FIXME: Support non-scalar expression
    vcode = code_formatter.visit(expr)
    conclusion = ['values[%d] = %s;' % (0, vcode)]

    code = [\
        '// Implementation of: %s' % str(expr),
        prelude, listing, conclusion,
        ]
    return code

# TODO: Support geometry from cell
#    /// Evaluate at given point in given cell 
#    virtual void eval(Array<double>& values, const Array<double>& x,
#                      const ufc::cell& cell) const;

# TODO: Support non-scalar expressions:
#    Expression();
#    Expression(uint dim);
#    Expression(uint dim0, uint dim1);
#    Expression(std::vector<uint> value_shape);

def compile_dolfin_expression(expr, name):

    assert expr.shape() == ()

    eval_body = compile_dolfin_expression_body(expr)

    eval_sig = 'virtual void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const'

    from ufl.algorithms import extract_coefficients
    coeffs = extract_coefficients(expr)
    #print coeffs # FIXME Use these to get object names

    coefficient_names = ['w%d' % k for k in range(len(coeffs))] # FIXME: get from object_names
    variables_code = ['boost::shared_ptr<dolfin::Function> %s;' % cname for cname in coefficient_names]

    class_body = [eval_sig, Block(eval_body)]
    if variables_code:
        class_body += ['', variables_code]

    code = [\
        'class Expression_%s: public dolfin::Expression' % name,
        '{', 'public:',
        Indented(class_body),
        '};',
        ]

    return code

test_header = """
#include <iostream>
#include <dolfin.h>
"""

test_footer = """
int main()
{
  uflacs::test::Expression_g g;
  return 0;
}
"""

def run_compile_dolfin(options, args):
    "Compile expressions from .ufl file(s) into dolfin C++ Expressions."
    from ufl.algorithms import load_ufl_file
    filenames = args
    namespace = 'uflacs'

    all_code = []

    # Generate code for each ufl file
    for fn in filenames:
        info("Loading file '%s'..." % (fn,))
        data = load_ufl_file(fn)
        file_code = []
        prefix, ext = os.path.splitext(os.path.basename(fn))
        if ext != '.ufl':
            print "Warning: expecting ufl file, got %s." % ext

        # Generate code for each expression in this file
        for k, expr in enumerate(data.expressions):
            name = data.object_names.get(id(expr), 'w%d' % k)
            expr_code = compile_dolfin_expression(expr, name)
            file_code.append(expr_code)
        
        # Wrap code from each file in its own namespace
        all_code.append(Namespace(prefix, file_code))

    # Wrap all code in the main namespace and render
    code = Namespace(namespace, all_code)
    code = format_code(code)

    # TODO: Write to file
    filecontents = '\n'.join((test_header, code, test_footer))
    print filecontents
    open('temp/main.cpp', 'w').write(filecontents)

    return 0
