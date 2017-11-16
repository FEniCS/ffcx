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
from ffc.quadrature.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])


def testProduct():
    "Test simple product instance."

    f_0 = format["float"](0)
    f_1 = format["float"](1)
    f0 = FloatValue(-2.0)
    f1 = FloatValue(3.0)
    f2 = FloatValue(0)
    f3 = FloatValue(-1)
    f4 = FloatValue(1)
    f5 = FloatValue(-0.5)
    f6 = FloatValue(2.0)
    s0 = Symbol("x", BASIS)
    s1 = Symbol("y", GEO)
    s2 = Symbol("z", GEO)

    p0 = Product([])
    p1 = Product([s0])
    p2 = Product([s0, s1])
    p3 = Product([f1, s0, s1])
    p4 = Product([s0, f2, s2])
    p5 = Product([s0, f0, s1, f1, s2])
    p6 = Product([s0, f3, s1])
    p7 = Product([s0, f4, s1]).expand().reduce_ops()
    p8 = Product([s0, f0, s2, f5])
    p9 = Product([s0, s1])
    p10 = Product([p0, p1])
    p11 = Product([f5, f0])
    p12 = Product([f6, f5])
    p13 = Product([f6, f5]).expand()
    p14 = Product([f1, f2])
    p_tmp = Product([f1])
    p_tmp.expand()
    p15 = Product([p_tmp, s0])

    assert repr(p0) == "Product([FloatValue(%s)])" % f_0
    assert repr(p1) == "Product([Symbol('x', BASIS)])"
    assert repr(p3) == "Product([FloatValue(%s), Symbol('x', BASIS), Symbol('y', GEO)])"\
        % format["float"](3)
    assert repr(p6) == "Product([FloatValue(-%s), Symbol('x', BASIS), Symbol('y', GEO)])" % f_1
    assert repr(p7) == "Product([Symbol('x', BASIS), Symbol('y', GEO)])"
    assert repr(p8) == "Product([Symbol('x', BASIS), Symbol('z', GEO)])"
    assert str(p2) == 'x*y'
    assert str(p4) == '%s' % f_0
    assert str(p5) == '-%s*x*y*z' % format["float"](6)
    assert str(p6) == ' - x*y'
    assert str(p7) == 'x*y'
    assert str(p8) == 'x*z'
    assert str(p9) == 'x*y'
    assert p0.val == 0
    assert str(p10) == '%s' % f_0
    assert str(p11) == '%s' % f_1
    assert str(p12) == '-%s' % f_1
    assert str(p13) == '-%s' % f_1
    assert repr(p14) == "Product([FloatValue(%s)])" % f_0
    assert repr(p14.expand()) == "FloatValue(%s)" % f_0

    assert p1 == p1
    assert p1 != p7
    assert p4 != p3
    assert p2 == p9
    assert p2 != p3

    assert p0.ops() == 0
    assert p1.ops() == 0
    assert p2.ops() == 1
    assert p3.ops() == 2
    assert p4.ops() == 0
    assert p5.ops() == 3
    assert p6.ops() == 1
    assert p7.ops() == 1
    assert p8.ops() == 1
    assert p9.ops() == 1
    assert p10.ops() == 0
    assert p14.ops() == 0

    # Test hash
    ll = [p3]
    d = {p3: 0}
    p10 = Product([f1, s0, s1])

    assert p3 in ll
    assert p3 in d
    assert p10 in ll
    assert p10 in d
    assert p2 not in ll
    assert p2 not in d
