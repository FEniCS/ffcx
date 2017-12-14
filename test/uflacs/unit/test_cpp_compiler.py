# -*- coding: utf-8 -*-
"""
Tests of generic C++ compilation code.
"""

# FIXME: These are all disabled, don't remember why, turn them on again!

import numpy


import ufl
from ufl import as_ufl
from ufl import *

import ffc.uflacs
# from ffc.uflacs.backends.toy.toy_compiler import compile_form


def format_variable(i):
    return "s[%d]" % i


def format_assignment(v, e):
    return "%s = %s;" % (v, e)


def format_code_lines(expressions, expression_dependencies, want_to_cache,
                      format_expression, format_variable, format_assignment):
    """FIXME: Test this.

    expressions             - array of length n with scalar valued expressions
    expression_dependencies - array of length n with
    want_to_cache           - array of length n with truth value for each expression whether we want to store it in a variable or not

    format_expression - a multifunction taking an expression as the first argument and code for each operand, returning code for the expression
    format_variable   - a function taking an int and formatting a variable name with it
    format_assignment - a function taking a lhs variable name and an rhs expression, returning an assignment statement.
    """
    def par(c):
        return "(%s)" % (c,)
    n = len(expressions)
    code_ref = numpy.empty(n, dtype=object)
    precedence = numpy.zeros(n, dtype=int)
    code_lines = []
    j = 0
    for i, e in enumerate(expressions):
        p = 0  # TODO: precedence of e
        deps = expression_dependencies[i]
        code_ops = [code_ref[d] for d in deps]
        code_ops = [par(c) if precedence[d] > p else c for c, d in zip(code_ops, deps)]
        code_e = format_expression(e, code_ops)
        if want_to_cache[i]:
            varname = format_variable(j)
            j += 1  # TODO: allocate free variable instead of just adding a new one
            assignment = format_assignment(varname, code_e)
            code_lines.append(assignment)
            code_ref[i] = varname
            precedence[i] = 9999  # TODO: highest # no () needed around variable
        else:
            code_ref[i] = code_e
            precedence[i] = p
        # TODO: deallocate variables somehow when they have been last used
    final_expressions = []  # TODO: Insert code expressions for given output expression indices here
    return code_lines, final_expressions

# class CppExpressionCompilerTest(UflTestCase):


def xtest_literal_zero_compilation():

    uexpr = as_ufl(0)

    lines, finals = compile_expression_lines(uexpr)

    expected_lines = []
    expected_finals = ['0']

    assert lines == expected_lines
    assert finals == expected_finals


def xtest_literal_int_compilation():

    uexpr = as_ufl(2)

    lines, finals = compile_expression_lines(uexpr)

    expected_lines = []
    expected_finals = ['2']

    assert lines == expected_lines
    assert finals == expected_finals


def xtest_literal_float_compilation():

    uexpr = as_ufl(2.56)

    lines, finals = compile_expression_lines(uexpr)

    expected_lines = []
    expected_finals = ['2.56']

    assert lines == expected_lines
    assert finals == expected_finals


def xtest_geometry_h_compilation():
    h = ufl.Circumradius(ufl.triangle)

    uexpr = h

    lines, finals = compile_expression_lines(uexpr)

    expected_lines = []
    expected_finals = ['h']

    assert lines == expected_lines
    assert finals == expected_finals


def xtest_terminal_sum_compilation():
    h = ufl.Circumradius(ufl.triangle)

    uexpr = h + 2.56

    lines, finals = compile_expression_lines(uexpr)

    expected_lines = []
    expected_finals = ['h + 2.56']

    assert lines == expected_lines
    assert finals == expected_finals


# class CppCompilerTest(UflTestCase):
import pytest


@pytest.fixture
def u():
    return Coefficient(FiniteElement("U", triangle, 1))


@pytest.fixture
def v():
    return Coefficient(VectorElement("U", triangle, 1))


@pytest.fixture
def w():
    return Coefficient(TensorElement("U", triangle, 1))


def test_fixtures(u, v, w):
    "Just checking that fixtures work!"
    assert u == Coefficient(FiniteElement("U", triangle, 1), count=u.count())
    assert v == Coefficient(VectorElement("U", triangle, 1), count=v.count())
    assert w == Coefficient(TensorElement("U", triangle, 1), count=w.count())


def xtest_cpp2_compile_scalar_literals():
    M = as_ufl(0) * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected

    M = as_ufl(3) * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected

    M = as_ufl(1.03) * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected


def xtest_cpp2_compile_geometry():
    M = CellVolume(triangle) * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected

    M = SpatialCoordinate(triangle)[0] * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected


def xtest_cpp2_compile_coefficients(u, v, w):
    M = u * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected

    M = v[0] * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected

    M = w[1, 0] * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected


def xtest_cpp2_compile_sums(u, v, w):
    M = (2 + u) * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected

    M = (v[1] + w[1, 1] + 3 + u) * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected


def xtest_cpp2_compile_products(u, v, w):
    M = (2 * u) * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected

    M = (v[1] * w[1, 1] * 3 * u) * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected


def xtest_cpp_compilation():
    M = u**2 / 2 * dx
    code = compile_form(M, 'unittest')
    print('\n', code)
    expected = 'TODO'
    assert code == expected
