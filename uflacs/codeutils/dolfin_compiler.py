
from ufl.classes import Terminal, Indexed, SpatialDerivative
from ufl.algorithms import Graph, preprocess_expression
from uflacs.codeutils.cpp_format import CppFormatterRules
from uflacs.codeutils.expr_formatter import ExprFormatter
from uflacs.codeutils.format_code import format_code, Block, Indented, Namespace


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
        basename = self.coefficient_names[o]
        # FIXME: Generate code for computing coefficient in prelude
        code = "%s[%d]" % (basename, o.count())
        if component:
            error("not implemented")
        if derivatives:
            error("not implemented")
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
    # Construct a specialized C++ expression formatter
    target_formatter = DolfinExpressionFormatter()
    cpp_formatter = CppFormatterRules(target_formatter)
    variables = {}
    expr_formatter = ExprFormatter(cpp_formatter, variables)

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
            vcode = expr_formatter.visit(v)
            listing.append('%s = %s; // %s' % (vname, vcode, str(v)))
            variables[v] = vname

    # Generate code prelude
    # FIXME: Compute coefficients
    # FIXME: Compute geometry
    prelude = ['double s[%d];' % k]

    # Generate assignments to values[] FIXME: Support non-scalar expression
    vcode = expr_formatter.visit(expr)
    conclusion = ['values[%d] = %s;' % (0, vcode)]

    code = [\
        '// Implementation of: %s' % str(expr),
        prelude, listing, conclusion,
        ]
    return code

def compile_dolfin_expression(expr, name, object_names):

    assert expr.shape() == ()

    eval_body = compile_dolfin_expression_body(expr)

    eval_sig = 'virtual void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const'

    from ufl.algorithms import extract_coefficients
    coeffs = extract_coefficients(expr)

    coefficient_names1 = [object_names.get(id(coeff), 'w%d' % coeff.count())\
                             for coeff in coeffs]
    coefficient_names2 = ['w%d' % k for k in range(len(coeffs))] # FIXME: get from object_names
    coefficient_names = coefficient_names1

    if 0: # debugging
        print
        print '\n'.join(map(str, zip(coefficient_names1, coeffs)))
        print
        print '\n'.join(map(str, zip(coefficient_names2, coeffs)))
        print

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
