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


def testSumOperators():
    "Test binary operators"

    f_0_5 = format["float"](0.5)
    f_1 = format["float"](1)
    f_2 = format["float"](2)
    f_3 = format["float"](3)
    f_6 = format["float"](6)
    f2 = FloatValue(2.0)
    fm3 = FloatValue(-3.0)

    x = Symbol("x", GEO)
    y = Symbol("y", GEO)
    z = Symbol("z", GEO)

    p0 = Product([f2, x])
    p1 = Product([x, y])

    S0 = Sum([x, y])
    S1 = Sum([x, z])

    F0 = Fraction(p0, y)

    # Test Sum '+'
    assert str(S0 + f2) == '(%s + x + y)' % f_2
    assert str(S0 + x) == '(x + x + y)'
    assert str(S0 + p0) == '(x + y + %s*x)' % f_2
    assert str(S0 + S0) == '(x + x + y + y)'
    assert str(S0 + F0) == '(x + y + %s*x/y)' % f_2

    # Test Sum '-'
    assert str(S0 - f2) == '(x + y-%s)' % f_2
    assert str(S0 - fm3) == '(x + y + %s)' % f_3
    assert str(S0 - x) == '(x + y - x)'
    assert str(S0 - p0) == '(x + y-%s*x)' % f_2
    assert str(S0 - Product([fm3, p0])) == '(x + y + %s*x)' % f_6
    assert str(S0 - S0) == '(x + y - (x + y))'
    assert str(S0 - F0) == '(x + y - %s*x/y)' % f_2

    # Test Sum '*'
    assert str(S0 * f2) == '(%s*x + %s*y)' % (f_2, f_2)
    assert str(S0 * x) == '(x*x + x*y)'
    assert str(S0 * p0) == '(%s*x*x + %s*x*y)' % (f_2, f_2)
    assert str(S0 * S0) == '(%s*x*y + x*x + y*y)' % f_2
    assert str(S0 * F0) == '(%s*x + %s*x*x/y)' % (f_2, f_2)

    # Test Sum '/'
    assert str(S0 / f2) == '(%s*x + %s*y)' % (f_0_5, f_0_5)
    assert str(S0 / x) == '(%s + y/x)' % f_1
    assert str(S0 / p0) == '(%s + %s*y/x)' % (f_0_5, f_0_5)
    assert str(S0 / p1) == '(%s/x + %s/y)' % (f_1, f_1)
    assert str(S0 / S0) == '(x + y)/(x + y)'
    assert str(S0 / S1) == '(x + y)/(x + z)'

    with pytest.raises(Exception):
        truediv(S0, FloatValue(0))
    with pytest.raises(Exception):
        truediv(S0, F0)
