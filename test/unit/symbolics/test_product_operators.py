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


def testProductOperators():
    "Test binary operators"

    f_0 = format["float"](0)
    f_2 = format["float"](2)
    f_4 = format["float"](4)

    f0 = FloatValue(0.0)
    f1 = FloatValue(1.0)
    f2 = FloatValue(2.0)
    fm1 = FloatValue(-1.0)
    fm3 = FloatValue(-3.0)

    x = Symbol("x", GEO)
    y = Symbol("y", GEO)
    z = Symbol("z", GEO)

    p0 = Product([f2, x])
    p1 = Product([x, y])
    p2 = Product([f2, z])
    p3 = Product([x, y, z])

    S0 = Sum([x, y])
    S1 = Sum([x, z])

    F0 = Fraction(f2, x)
    F1 = Fraction(x, y)
    F2 = Fraction(x, S0)
    F3 = Fraction(x, y)
    F4 = Fraction(p0, y)

    # Test Product '+'
    assert str(p0 + f2) == '(%s + %s*x)' % (f_2, f_2)
    assert str(p0 + x) == '%s*x' % format["float"](3)
    assert str(p0 + y) == '(y + %s*x)' % f_2
    assert str(p0 + p0) == '%s*x' % f_4
    assert str(p0 + p1) == '(%s*x + x*y)' % f_2
    assert p0 + Product([fm1, x]) == x
    assert Product([fm1, x]) + x == f0
    assert str(x + Product([fm1, x])) == '%s' % f_0
    assert str(p0 + S0) == '(x + y + %s*x)' % f_2
    assert str(p0 + F3) == '(%s*x + x/y)' % f_2

    # Test Product '-'
    assert str(p0 - f2) == '(%s*x-%s)' % (f_2, f_2)
    assert str(p0 - x) == 'x'
    assert str(p0 - y) == '(%s*x - y)' % f_2
    assert str(p0 - p0) == '%s' % f_0
    assert str(p0 - p1) == '(%s*x - x*y)' % f_2
    assert str(p0 - S0) == '(%s*x - (x + y))' % f_2
    assert str(p0 - F3) == '(%s*x - x/y)' % f_2

    # Test Product '*', only need to test float, symbol and product.
    # Sum and fraction are handled by 'other'
    assert str(p0 * f0) == '%s' % f_0
    assert str(p0 * fm3) == '-%s*x' % format["float"](6)
    assert str(p0 * y) == '%s*x*y' % f_2
    assert str(p0 * p1) == '%s*x*x*y' % f_2

    # Test Product '/'
    assert str(Product([f0, x]) / x) == '%s' % f_0
    assert str(p0 / S0) == '%s*x/(x + y)' % f_2
    assert p1 / y == x
    assert p1 / p2 == Fraction(Product([p1, FloatValue(0.5)]), z)
    assert p1 / z == Fraction(p1, z)
    assert p0 / Product([f2, p1]) == Fraction(f1, y)
    assert p1 / p0 == Product([FloatValue(0.5), y])
    assert p1 / p1 == f1
    assert p1 / p3 == Fraction(f1, z)
    assert str(p1 / p3) == '%s/z' % format["float"](1)

    with pytest.raises(Exception):
        truediv(p0, f0)
    with pytest.raises(Exception):
        truediv(p0, F0)
