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


def testFractionOperators():
    "Test binary operators"

    f_0 = format["float"](0)
    f_1 = format["float"](1)
    f_2 = format["float"](2)
    f_5 = format["float"](5)

    f2 = FloatValue(2.0)
    fm3 = FloatValue(-3.0)

    x = Symbol("x", GEO)
    y = Symbol("y", GEO)

    p0 = Product([f2, x])
    p1 = Product([x, y])

    S0 = Sum([x, y])

    F0 = Fraction(f2, y)
    F1 = Fraction(x, y)
    F2 = Fraction(x, S0)
    F3 = Fraction(x, y)
    F4 = Fraction(p0, y)
    F5 = Fraction(Product([fm3, x]), y)

    # Test Fraction '+'
    assert str(F0 + f2) == '(%s + %s/y)' % (f_2, f_2)
    assert str(F1 + x) == '(x + x/y)'
    assert str(F1 + p0) == '(%s*x + x/y)' % f_2
    assert str(F1 + S0) == '(x + y + x/y)'
    assert str(F1 + F3) == '%s*x/y' % f_2
    assert str(F0 + F1) == '(%s + x)/y' % f_2
    assert str(F2 + F4) == '(%s*x/y + x/(x + y))' % f_2

    # Test Fraction '-'
    assert str(F0 - f2) == '(%s/y-%s)' % (f_2, f_2)
    assert str(F1 - x) == '(x/y - x)'
    assert str(F1 - p0) == '(x/y-%s*x)' % f_2
    assert str(F1 - S0) == '(x/y - (x + y))'
    assert str(F1 - F3) == '%s' % f_0
    assert str(F4 - F1) == 'x/y'
    assert str(F4 - F5) == '%s*x/y' % f_5
    assert str(F0 - F1) == '(%s - x)/y' % f_2
    assert str(F2 - F4) == '(x/(x + y) - %s*x/y)' % f_2

    # Test Fraction '*'
    assert str(F1 * f2) == '%s*x/y' % f_2
    assert str(F1 * x) == 'x*x/y'
    assert str(F1 * p1) == 'x*x'
    assert str(F1 * S0) == '(x + x*x/y)'
    assert repr(F1 * S0) == repr(Sum([x, Fraction(Product([x, x]), y)]))
    assert str(F1 * F0) == '%s*x/(y*y)' % f_2

    # Test Fraction '/'
    assert str(F0 / f2) == '%s/y' % f_1
    assert str(F1 / x) == '%s/y' % f_1
    assert str(F4 / p1) == '%s/(y*y)' % f_2
    assert str(F4 / x) == '%s/y' % f_2
    assert str(F2 / y) == 'x/(x*y + y*y)'
    assert str(F0 / S0) == '%s/(x*y + y*y)' % f_2

    with pytest.raises(Exception):
        truediv(F0 / F0)
