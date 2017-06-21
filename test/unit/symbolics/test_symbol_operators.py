# -*- coding: utf-8 -*-
# Copyright (C) 2010 Kristian B. Oelgaard
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.

import pytest
import time

# FFC modules
from ffc.quadrature.reduce_operations import operation_count, expand_operations, reduce_operations
from ffc.quadrature.symbolics import *
from ffc.quadrature.sumobj import _group_fractions
from ffc.quadrature.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])
from ffc.log import error, push_level, pop_level, CRITICAL


def testSymbolOperators():
    "Test binary operators"

    f_0 = format["float"](0)
    f_1 = format["float"](1)
    f_2 = format["float"](2)
    f_3 = format["float"](3)
    f_0_5 = format["float"](0.5)
    f0 = FloatValue(0.0)
    f2 = FloatValue(2.0)
    fm1 = FloatValue(-1.0)
    fm3 = FloatValue(-3.0)

    x = Symbol("x", GEO)
    y = Symbol("y", GEO)
    z = Symbol("z", GEO)

    p0 = Product([f2, x])
    p1 = Product([x, y])
    p2 = Product([f2, z])
    p3 = Product([y, x, z])

    S0 = Sum([x, y])
    S1 = Sum([x, z])

    F0 = Fraction(f2, y)
    F1 = Fraction(x, y)
    F2 = Fraction(x, S0)
    F3 = Fraction(x, y)
    F4 = Fraction(p0, y)
    F5 = Fraction(fm3, y)

    # Test Symbol '+'
    assert str(x + f2) == '(%s + x)' % f_2
    assert str(x + x) == '%s*x' % f_2
    assert str(x + y) == '(x + y)'
    assert str(x + p0) == '%s*x' % f_3
    assert str(x + p1) == '(x + x*y)'
    assert str(x + S0) == '(x + x + y)'
    assert str(x + F0) == '(x + %s/y)' % f_2

    # Test Symbol '-'
    assert str(x - f2) == '(x-%s)' % f_2
    assert str(x - x) == '%s' % f_0
    assert str(x - y) == '(x - y)'
    assert str(x - p0) == ' - x'
    assert str(x - p1) == '(x - x*y)'
    assert str(x - S0) == '(x - (x + y))'
    assert str(x - F5) == '(x - -%s/y)' % f_3

    # Test Symbol '*', only need to test float, symbol and product. Sum and
    # fraction are handled by 'other'
    assert str(x * f2) == '%s*x' % f_2
    assert str(x * y) == 'x*y'
    assert str(x * p1) == 'x*x*y'

    # Test Symbol '/'
    assert str(x / f2) == '%s*x' % f_0_5
    assert str(x / x) == '%s' % f_1
    assert str(x / y) == 'x/y'
    assert str(x / S0) == 'x/(x + y)'
    assert str(x / p0) == '%s' % f_0_5
    assert str(y / p1) == '%s/x' % f_1
    assert str(z / p0) == '%s*z/x' % f_0_5
    assert str(z / p1) == 'z/(x*y)'
    with pytest.raises(Exception):
        truediv(x, F0)
    with pytest.raises(Exception):
        truediv(y, FloatValue(0))
