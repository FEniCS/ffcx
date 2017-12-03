# -*- coding: utf-8 -*-
"""
Tests of dolfin Expression formatting.
"""

import pytest

import ufl
from ufl.constantvalue import as_ufl

import numpy

import hashlib

from ffc.uflacs.params import default_parameters
from ffc.uflacs.backends.dolfin.expression import format_dolfin_expression
from ffc.uflacs.backends.dolfin.dolfin_compiler import compile_dolfin_expression_body


# FIXME: Make tests with dolfin optional
@pytest.fixture
def dolfin():
    import dolfin
    return dolfin


# For test comparisons
dolfin_expression_minimal = """\
class MyExpression: public Expression
{
public:

    virtual void eval(Array<double>& values,
                      const Array<double>& x) const
    {

    }
};"""

dolfin_expression_with_constructor1 = """\
class MyExpression: public Expression
{
public:
    MyExpression(): Expression(2) {}

    virtual void eval(Array<double>& values,
                      const Array<double>& x) const
    {

    }
};"""

dolfin_expression_with_constructor2 = """\
class MyExpression: public Expression
{
public:
    MyExpression(): Expression(2, 3) {}

    virtual void eval(Array<double>& values,
                      const Array<double>& x) const
    {

    }
};"""

dolfin_expression_with_members = """\
class MyExpression: public Expression
{
public:

    virtual void eval(Array<double>& values,
                      const Array<double>& x) const
    {

    }

    double c;
    boost::shared_ptr<MeshFunction<double> > mf;
    boost::shared_ptr<GenericFunction> gf;
    boost::shared_ptr<Function> f;
};"""

dolfin_expression_with_members2 = """\
class MyExpression: public Expression
{
public:

    virtual void eval(Array<double>& values,
                      const Array<double>& x) const
    {

    }

    double c0;
    double c1;
    boost::shared_ptr<MeshFunction<double> > mf0;
    boost::shared_ptr<MeshFunction<double> > mf1;
    boost::shared_ptr<GenericFunction> gf0;
    boost::shared_ptr<GenericFunction> gf1;
    boost::shared_ptr<Function> f0;
    boost::shared_ptr<Function> f1;
};"""

dolfin_expression_with_cell = """\
class MyExpression: public Expression
{
public:

    void eval(Array<double>& values, const Array<double>& x,
              const ufc::cell& ufc_cell) const
    {
        // Emulate tabulate_tensor arguments
        const ufc::cell& c = ufc_cell;
        Array<double> & A = values;

        //dolfin_assert(ufc_cell.local_facet >= 0);

        values[0] = 1.23;
    }
};"""


other_stuff = """
        // Getting cell data
        const uint gd = cell.geometric_dimension;
        const uint td = cell.topological_dimension;
        const uint cell_index = cell.entity_indices[td][0];

        // Tabulation of mesh function values
        double mf0v = (*mf0)[cell_index];

        // Evaluation of functions
        Array<double> w0v(1);
        w0->eval(w0v, x, cell);

        // Writing output values
        values[0] = w0v[0]*mf0v;
        values[1] = w0v[0]+mf0v;
"""

#class DolfinExpressionFormatterTest(UflTestCase):
#"""Tests of Dolfin C++ Expression formatting utilities."""

def test_dolfin_expression_with_cell(dolfin):
    e = dolfin.Expression(cppcode=dolfin_expression_with_cell)
    assert e.ufl_shape == ()
    values = numpy.zeros((1,))
    x = numpy.asarray((0.2, 0.2))
    mesh = dolfin.UnitSquareMesh(1, 1)
    cell = dolfin.Cell(mesh, 0)
    e.eval_cell(values, x, cell)
    assert values[0] == 1.23
    #assert e((1.0,2.0)) == 1.23 # FAILS! Need dolfin to dispatch this to eval_cell with the proper Cell found from x for this to work.

def test_dolfin_expression_constructors():
    code = format_dolfin_expression(classname="MyExpression",
                                    shape=())
    assert code == dolfin_expression_minimal
    code = format_dolfin_expression(classname="MyExpression",
                                    shape=(2,))
    assert code == dolfin_expression_with_constructor1
    code = format_dolfin_expression(classname="MyExpression",
                                    shape=(2, 3))
    assert code == dolfin_expression_with_constructor2

def test_dolfin_expression_members():
    code = format_dolfin_expression(classname="MyExpression",
                                    shape=(),
                                    constants=('c',),
                                    mesh_functions=('mf',),
                                    generic_functions=('gf',),
                                    functions=('f',),
                                    )
    assert code == dolfin_expression_with_members
    code = format_dolfin_expression(classname="MyExpression",
                                    shape=(),
                                    constants=('c0', 'c1'),
                                    mesh_functions=('mf0', 'mf1'),
                                    generic_functions=('gf0', 'gf1'),
                                    functions=('f0', 'f1'),
                                    )
    assert code == dolfin_expression_with_members2

def test_explicit_dolfin_expression_compilation(dolfin):
    code = format_dolfin_expression(classname="MyExpression",
                                    shape=(2,),
                                    eval_body=['values[0] = x[0];',
                                               'values[1] = x[1];'],
                                    constants=('c0', 'c1'),
                                    mesh_functions=('mf0', 'mf1'),
                                    generic_functions=('gf0', 'gf1'),
                                    functions=('f0', 'f1'),
                                    )
    expr = dolfin.Expression(cppcode=code)
    assert hasattr(expr, 'c0')
    assert hasattr(expr, 'c1')
    assert hasattr(expr, 'mf0')
    assert hasattr(expr, 'mf1')
    assert hasattr(expr, 'gf0')
    assert hasattr(expr, 'gf1')
    assert hasattr(expr, 'f0')
    assert hasattr(expr, 'f1')
    #dolfin.plot(expr, mesh=dolfin.UnitSquareMesh(10,10), interactive=True)


def flattened_nonempty_lines(lines):
    res = []
    for l in lines:
        if isinstance(l, list):
            res.extend(flattened_nonempty_lines(l))
        elif l:
            res.append(l)
        else:
            pass
    return res

#class Ufl2DolfinExpressionCompilerTest(UflTestCase):
#"""Tests of Dolfin C++ Expression compilation from UFL expressions."""

def check_dolfin_expression_compilation(uexpr, expected_lines, expected_values, members={}):
    parameters = default_parameters()
    # Compile expression
    compiled_lines, member_names = compile_dolfin_expression_body(uexpr, parameters)

    # Check expected compilation output
    if flattened_nonempty_lines(compiled_lines) != expected_lines:
        print('\n'*5)
        print("Upcoming failure, expected lines:")
        print('\n'.join(expected_lines))
        print("\nActual lines:")
        print('\n'.join(flattened_nonempty_lines(compiled_lines)))
        print('\n'*5)
    assert flattened_nonempty_lines(compiled_lines) == expected_lines

    # Add debugging lines to print scalar expressions at end of Expression eval call
    if 0:
        printme = ["x[0]", "values[0]"]
        compiled_lines += ['std::cout '] + ['       << "%s = " << %s << std::endl' % (pm, pm) for pm in printme] + ['<< std::endl;']

    # Wrap it in a dolfin::Expression class
    shape = uexpr.ufl_shape
    name = "MyExpression_%s" % hashlib.md5(str(hash(uexpr))).hexdigest()
    code = format_dolfin_expression(classname=name,
                                    shape=shape,
                                    eval_body=compiled_lines,
                                    **member_names)
    if 0:
        print("SKIPPING REST OF TEST")
        return

    # Try to compile it with DOLFIN
    import dolfin
    dexpr = dolfin.Expression(cppcode=code)

    # Connect compiled dolfin::Expression object with dolfin::Function instances
    for name, value in members.items():
        setattr(dexpr, name, value)

    # Evaluate and assert compiled value!
    for x, v in expected_values:
        u = dexpr(x) # TODO: Not sure if this works
        if not hasattr(u, '__len__'):
            u = (u,)
        diff = numpy.array(u) - numpy.array(v)
        if numpy.dot(diff, diff) > 1e-13:
            print("u =", u)
            print("v =", v)
            print("diff =", diff)
        assert numpy.dot(diff, diff) < 1e-14

def test_dolfin_expression_compilation_of_scalar_literal():
    # Define some literal ufl expression
    uexpr = as_ufl(3.14)

    # Define expected output from compilation
    expected_lines = ['double s[1];',
                      's[0] = 3.14;',
                      'values[0] = s[0];']

    # Define expected evaluation values: [(x,value), (x,value), ...]
    expected_values = [((0.0, 0.0), (3.14,)),
                       ((0.6, 0.7), (3.14,)),
                       ]

    # Execute all tests
    check_dolfin_expression_compilation(uexpr, expected_lines, expected_values)

def test_dolfin_expression_compilation_of_vector_literal():
    # Define some literal ufl expression
    uexpr = ufl.as_vector((1.23, 7.89))

    # Define expected output from compilation
    expected_lines = ['double s[2];',
                      's[0] = 1.23;',
                      's[1] = 7.89;',
                      'values[0] = s[0];',
                      'values[1] = s[1];']

    # Define expected evaluation values: [(x,value), (x,value), ...]
    expected_values = [((0.0, 0.0), (1.23, 7.89)),
                       ((0.6, 0.7), (1.23, 7.89)),
                       ]

    # Execute all tests
    check_dolfin_expression_compilation(uexpr, expected_lines, expected_values)

def test_dolfin_expression_compilation_of_x():
    # Define some ufl expression
    x = ufl.SpatialCoordinate(ufl.triangle)
    uexpr = 2*x

    # Define expected output from compilation
    expected_lines = ['double s[2];',
                      's[0] = 2 * x[0];',
                      's[1] = 2 * x[1];',
                      'values[0] = s[0];',
                      'values[1] = s[1];']

    # Define expected evaluation values: [(x,value), (x,value), ...]
    expected_values = [((0.0, 0.0), (0.0, 0.0)),
                       ((0.6, 0.7), (1.2, 1.4)),
                       ]

    # Execute all tests
    check_dolfin_expression_compilation(uexpr, expected_lines, expected_values)

def test_dolfin_expression_compilation_of_coefficient(dolfin):

    # Define some ufl expression with PyDOLFIN coefficients
    mesh = dolfin.UnitSquareMesh(3, 3)
    # Using quadratic element deliberately for accuracy
    V = dolfin.FunctionSpace(mesh, "CG", 2)
    u = dolfin.Function(V)
    u.interpolate(dolfin.Expression("x[0]*x[1]"))
    uexpr = 2*u

    # Define expected output from compilation
    expected_lines = ['double s[1];',
                      'Array<double> v_w0(1);',
                      'w0->eval(v_w0, x);',
                      's[0] = 2 * v_w0[0];',
                      'values[0] = s[0];']

    # Define expected evaluation values: [(x,value), (x,value), ...]
    expected_values = [((0.0, 0.0), (0.0,)),
                       ((0.6, 0.7), (2*0.6*0.7,)),
                       ]

    # Execute all tests
    check_dolfin_expression_compilation(uexpr, expected_lines, expected_values, members={'w0':u})

def test_dolfin_expression_compilation_of_algebraic_operators(dolfin):

    # Define some PyDOLFIN coefficients
    mesh = dolfin.UnitSquareMesh(3, 3)
    # Using quadratic element deliberately for accuracy
    V = dolfin.FunctionSpace(mesh, "CG", 2)
    u = dolfin.Function(V)
    u.interpolate(dolfin.Expression("x[0]*x[1]"))

    # Define ufl expression with algebraic operators
    uexpr = (2*u + 3.14*u**2) / 4

    # Define expected output from compilation
    # TODO: Test with and without variables
    # Fully split version:
    expected_lines = ['double s[4];',
                      'Array<double> v_w0(1);',
                      'w0->eval(v_w0, x);',
                      's[0] = 2 * v_w0[0];',
                      's[1] = pow(v_w0[0], 2);',
                      's[2] = 3.14 * s[1];',
                      's[3] = s[2] / 4;',
                      'values[0] = s[3];']
    # Oneliner version (no reuse in this expression):
    expected_lines = ['double s[1];',
                      'Array<double> v_w0(1);',
                      'w0->eval(v_w0, x);',
                      's[0] = (3.14 * pow(v_w0[0], 2) + 2 * v_w0[0]) / 4;',
                      'values[0] = s[0];']

    # Define expected evaluation values: [(x,value), (x,value), ...]
    x, y = 0.6, 0.7; expected = (2*x*y + 3.14*(x*y)**2) / 4
    expected_values = [((0.0, 0.0), (0.0,)),
                       ((x, y), (expected,)),
                       ]

    # Execute all tests
    check_dolfin_expression_compilation(uexpr, expected_lines, expected_values, members={'w0':u})


def test_dolfin_expression_compilation_of_math_functions(dolfin):

    # Define some PyDOLFIN coefficients
    mesh = dolfin.UnitSquareMesh(3, 3)
    # Using quadratic element deliberately for accuracy
    V = dolfin.FunctionSpace(mesh, "CG", 2)
    u = dolfin.Function(V)
    u.interpolate(dolfin.Expression("x[0]*x[1]"))
    w0 = u

    # Define ufl expression with math functions
    v = abs(ufl.cos(u))/2 + 0.02
    uexpr = ufl.sin(u) + ufl.tan(v) + ufl.exp(u) + ufl.ln(v) + ufl.atan(v) + ufl.acos(v) + ufl.asin(v)

    #print dolfin.assemble(uexpr**2*dolfin.dx, mesh=mesh) # 11.7846508409

    # Define expected output from compilation
    ucode = 'v_w0[0]'
    vcode = '0.02 + fabs(cos(v_w0[0])) / 2'
    funcs = 'asin(%(v)s) + (acos(%(v)s) + (atan(%(v)s) + (log(%(v)s) + (exp(%(u)s) + (sin(%(u)s) + tan(%(v)s))))))'
    oneliner = funcs % {'u':ucode, 'v':vcode}

    # Oneliner version (ignoring reuse):
    expected_lines = ['double s[1];',
                      'Array<double> v_w0(1);',
                      'w0->eval(v_w0, x);',
                      's[0] = %s;' % oneliner,
                      'values[0] = s[0];']
    #cppcode = format_dolfin_expression(classname="DebugExpression", shape=(), eval_body=expected_lines)
    #print '-'*100
    #print cppcode
    #print '-'*100
    #dolfin.plot(dolfin.Expression(cppcode=cppcode, mesh=mesh))
    #dolfin.interactive()

    # Split version (handles reuse of v, no other reuse):
    expected_lines = ['double s[2];',
                      'Array<double> v_w0(1);',
                      'w0->eval(v_w0, x);',
                      's[0] = %s;' % (vcode,),
                      's[1] = %s;' % (funcs % {'u':ucode,'v':'s[0]'},),
                      'values[0] = s[1];']

    # Define expected evaluation values: [(x,value), (x,value), ...]
    import math
    x, y = 0.6, 0.7
    u = x*y
    v = abs(math.cos(u))/2 + 0.02
    v0 = .52
    expected0 = math.tan(v0) + 1 + math.log(v0) + math.atan(v0) + math.acos(v0) + math.asin(v0)
    expected = math.sin(u) + math.tan(v) + math.exp(u) + math.log(v) + math.atan(v) + math.acos(v) + math.asin(v)
    expected_values = [((0.0, 0.0), (expected0,)),
                       ((x, y), (expected,)),
                       ]

    # Execute all tests
    check_dolfin_expression_compilation(uexpr, expected_lines, expected_values, members={'w0':w0})

def test_dolfin_expression_compilation_of_dot_with_index_notation(dolfin):

    # Define some PyDOLFIN coefficients
    mesh = dolfin.UnitSquareMesh(3, 3)
    # Using quadratic element deliberately for accuracy
    V = dolfin.VectorFunctionSpace(mesh, "CG", 2)
    u = dolfin.Function(V)
    u.interpolate(dolfin.Expression(("x[0]", "x[1]")))

    # Define ufl expression with the simplest possible index notation
    uexpr = ufl.dot(u, u) # Needs expand_compounds, not testing here

    uexpr = u[0]*u[0] + u[1]*u[1] # Works

    i, = ufl.indices(1)
    uexpr = u[i]*u[i] # Works

    # Define expected output from compilation
    # Oneliner version (no reuse in this expression):
    xc, yc = ('v_w0[0]', 'v_w0[1]')
    expected_lines = ['double s[1];',
                      'Array<double> v_w0(2);',
                      'w0->eval(v_w0, x);',
                      's[0] = pow(%(x)s, 2) + pow(%(y)s, 2);' % {'x':xc, 'y':yc},
                      'values[0] = s[0];']

    # Define expected evaluation values: [(x,value), (x,value), ...]
    x, y = (0.6, 0.7)
    expected = (x*x+y*y)
    expected_values = [((0.0, 0.0), (0.0,)),
                       ((x, y), (expected,)),
                       ]

    # Execute all tests
    check_dolfin_expression_compilation(uexpr, expected_lines, expected_values, members={'w0':u})

def test_dolfin_expression_compilation_of_index_notation(dolfin):

    # Define some PyDOLFIN coefficients
    mesh = dolfin.UnitSquareMesh(3, 3)
    # Using quadratic element deliberately for accuracy
    V = dolfin.VectorFunctionSpace(mesh, "CG", 2)
    u = dolfin.Function(V)
    u.interpolate(dolfin.Expression(("x[0]", "x[1]")))

    # Define ufl expression with index notation
    i, j = ufl.indices(2)
    uexpr = ((u[0]*u + u[1]*u)[i]*u[j])*(u[i]*u[j])

    # Define expected output from compilation
    # Split version (reusing all product terms correctly!):
    xc, yc = ('v_w0[0]', 'v_w0[1]')
    svars = {'a':'s[2]', 'b':'s[4]', 'xx':'s[0]', 'yy':'s[3]', 'xy':'s[1]', 'x':xc, 'y':yc}
    expected_lines = ['double s[6];',
                      'Array<double> v_w0(2);',
                      'w0->eval(v_w0, x);',
                      's[0] = pow(%(x)s, 2);' % svars,
                      's[1] = %(x)s * %(y)s;' % svars,
                      's[2] = %(xx)s + %(xy)s;' % svars,
                      's[3] = pow(%(y)s, 2);' % svars,
                      's[4] = %(yy)s + %(xy)s;' % svars,
                      's[5] = (%(xx)s * (%(x)s * %(a)s) + %(xy)s * (%(x)s * %(b)s)) + (%(yy)s * (%(y)s * %(b)s) + %(xy)s * (%(y)s * %(a)s));' % svars,
                      'values[0] = s[5];']
    # Making the above line correct was a nightmare because of operator sorting in ufl... Hoping it stays stable!

    # Define expected evaluation values: [(x,value), (x,value), ...]
    x, y = (0.6, 0.7)
    a = x*x+y*x
    b = x*y+y*y
    expected = (a*x)*(x*x) + (a*y)*(x*y) + (b*x)*(y*x) + (b*y)*(y*y)
    expected_values = [((0.0, 0.0), (0.0,)),
                       ((x, y), (expected,)),
                       ]

    # Execute all tests
    check_dolfin_expression_compilation(uexpr, expected_lines, expected_values, members={'w0':u})
