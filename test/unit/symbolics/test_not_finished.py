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
from ffc.quadrature.symbolics import *
from ffc.quadrature.sumobj import _group_fractions
from ffc.quadrature.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])


def testNotFinished():
    "Stuff that would be nice to implement."

    f_1 = format["float"](1)
    f_2 = format["float"](2)
    f_4 = format["float"](4)
    f_8 = format["float"](8)

    f0 = FloatValue(4)
    f1 = FloatValue(2)
    f2 = FloatValue(8)
    s0 = Symbol("x", GEO)
    s1 = Symbol("y", GEO)
    s2 = Symbol("z", GEO)
    a = Symbol("a", GEO)
    b = Symbol("b", GEO)
    c = Symbol("c", GEO)

    # Aux. expressions
    p0 = Product([f1, s0])
    p1 = Product([f2, s1])
    p2 = Product([s0, s1])

    F0 = Fraction(f0, s0)

    S0 = Sum([p0, p1])
    S1 = Sum([s0, p2])
    S2 = Sum([FloatValue(1), s1])
    S3 = Sum([F0, F0])

    # Thing to be implemented
    e0 = f0 / S0
    e1 = s0 / S1
    e2 = S2 / S1
    e3 = _group_fractions(S3)
    e4 = Sum([Fraction(f1 * s0, a * b * c), Fraction(s0, a * b)]).expand().reduce_ops()

    # Tests that pass the current implementation
    assert str(e0) == '%s/(%s*x + %s*y)' % (f_4, f_2, f_8)
    assert str(e1) == 'x/(x + x*y)'
    assert str(e2) == '(%s + y)/(x + x*y)' % f_1
    assert str(e3) == '%s/x' % f_8
    assert str(e4) == 'x*(%s/(a*b) + %s/(a*b*c))' % (f_1, f_2)

    # Tests that should pass in future implementations (change NotEqual to Equal)
    assert str(e0) != '%s/(x + %s*y)' % (f_2, f_4)
    assert str(e1) != '%s/(%s + y)' % (f_1, f_1)
    assert str(e2) != '%s/x' % f_1
    assert str(e4) != 'x*(%s/c + %s)/(a*b)' % (f_2, f_1)

    # TODO: Would it be a good idea to reduce expressions
    # wrt. var_type without first expanding?
    E0 = Product([Sum([Product([Symbol('B0', BASIS), Product([Symbol('B1', BASIS), Sum([s0]), Sum([s0])])]),
                       Product([Symbol('B0', BASIS), Symbol('B1', BASIS)])])])
    Er0 = E0.reduce_vartype(BASIS)
    Ex0 = E0.expand().reduce_vartype(BASIS)
    assert Ex0[0][1] != Er0[0][1].expand()

    # Both of these reductions should work at the same time
    # 1) 2/(x/(a+b) + y/(a+b)) --> 2(a+b)/(x+y)
    # 2) 2/(x + y/(a+b)) --> no reduction, or if divisions are more expensive
    # 3) 2/(x + y/(a+b)) --> 2(a+b)/((a+b)x + y)
