
from ufl.classes import Terminal, Indexed, SpatialDerivative
from ufl.algorithms import Graph, preprocess_expression
from uflacs.codeutils.cpp_format import CppFormatterRules
from uflacs.codeutils.expr_formatter import ExprFormatter
from uflacs.codeutils.format_code import format_code, Block, Indented, Namespace, Class


# TODO: Support geometry from cell
#    /// Evaluate at given point in given cell 
#    virtual void eval(Array<double>& values, const Array<double>& x,
#                      const ufc::cell& cell) const;

# TODO: Support non-scalar expressions:
#    Expression();
#    Expression(uint dim);
#    Expression(uint dim0, uint dim1);
#    Expression(std::vector<uint> value_shape);

# TODO: Generate a test code for validation of Expression?
test_template = """
#include <iostream>
#include <cmath>
#include <dolfin.h>
#include <{path}.h>

int main()
{
  uflacs::{namespace}::{classname} expr;
  // ...
  return 0;
}
"""


class DolfinExpressionFormatter(object):
    """Cpp formatter class for dolfin::Expression generation."""

    def __init__(self, expr_data): # TODO: Take in what's needed instead of expr_data
        self.coefficient_names = expr_data.coefficient_names

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

    def coefficient(self, o, component=(), derivatives=()):
        k = o.count()
        basename = self.coefficient_names[k]

        code = "%s[%d]" % (basename, k)
        if component:
            error("not implemented")
        if derivatives:
            error("not implemented")

        # FIXME: Register that this value is required,
        # then generate code for computing coefficient in prelude

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

def compile_dolfin_expression_body(expr_data):
    # Construct a specialized C++ expression formatter
    target_formatter = DolfinExpressionFormatter(expr_data)
    cpp_formatter = CppFormatterRules(target_formatter)
    variables = {}
    expr_formatter = ExprFormatter(cpp_formatter, variables)

    # Get preprocessed expression and build computational graph
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
            vcode = expr_formatter.visit(v)
            listing.append('%s = %s; // %s' % (vname, vcode, str(v)))
            variables[v] = vname

    # Generate code prelude
    # FIXME: Compute coefficients
    # FIXME: Compute geometry
    prelude = ['double s[%d];' % k] if k else []

    # Generate assignments to values[] FIXME: Support non-scalar expression
    vcode = expr_formatter.visit(expr)
    assignment = ['values[%d] = %s;' % (0, vcode)]

    code = [\
        '// Implementation of: %s' % str(expr),
        prelude, listing, assignment,
        ]
    return code

def compile_dolfin_expression(expr, name, object_names):

    assert expr.shape() == () # TODO: support this

    expr_data = preprocess_expression(expr, object_names=object_names)

    classname = 'Expression_%s' % name

    # Create a member for each coefficient
    variables_code = ['boost::shared_ptr<dolfin::Function> %s;' % cname for cname in expr_data.coefficient_names]

    # Choose constructor based on value shape
    r = len(expr.shape())
    if r == 0:
        constructor = '' # default constructor is fine
    elif r == 1:
        constructor = (classname, '(dolfin::uint dim): dolfin::Expression(dim) {}')
    elif r == 2:
        constructor = (classname, '(dolfin::uint dim0, dolfin::uint dim1): ',
                       'dolfin::Expression(dim0, dim1) {}')
    else:
        constructor = (classname, '(std::vector<dolfin::uint> value_shape): ',
                       'dolfin::Expression(value_shape) {}')

    eval_sig = 'virtual void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const'

    eval_body = compile_dolfin_expression_body(expr_data)

    class_body = [constructor, eval_sig, Block(eval_body)]
    if variables_code:
        class_body += ['', variables_code]

    code = Class(classname, 'dolfin::Expression', public_body=class_body)
    code = format_code(code)
    return code, classname


# TODO: Move this to dolfin compiler file
def compile_dolfin_expressions_header(data, prefix):
    from uflacs.codeutils.dolfin_compiler import compile_dolfin_expression

    includes = ['#include <iostream>', '#include <cmath>', '#include <dolfin.h>']

    # Generate code for each expression in this file
    file_code = []
    for k, expr in enumerate(data.expressions):
        name = data.object_names.get(id(expr), 'w%d' % k)
        expr_code, expr_classname = compile_dolfin_expression(expr, name, data.object_names)
        file_code.append(expr_code)

    # Wrap code from each file in its own namespace
    define = 'UFLACS_' + prefix + '_INCLUDED'
    preguard = [('#ifndef ', define), ('#define ', define)]
    postguard = '#endif'
    code = format_code([preguard, '', includes, '', Namespace(prefix, file_code), '', postguard])
    return code
