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


def testSum():
    "Test simple sum instance."

    f_0 = format["float"](0)
    f_1 = format["float"](1)
    f_2 = format["float"](2)
    f_3 = format["float"](3)

    f0 = FloatValue(-2.0)
    f1 = FloatValue(3.0)
    f2 = FloatValue(0)
    s0 = Symbol("x", BASIS)
    s1 = Symbol("y", GEO)
    s2 = Symbol("z", GEO)

    S0 = Sum([])
    S1 = Sum([s0])
    S2 = Sum([s0, s1])
    S3 = Sum([s0, s0])
    S4 = Sum([f0, s0])
    S5 = Sum([s0, f0, s0])
    S6 = Sum([s0, f0, s0, f1])
    S7 = Sum([s0, f0, s1, f2])
    S8 = Sum([s0, f1, s0])
    S9 = Sum([f0, f0, f0, f1, f1, s1])
    S10 = Sum([s1, s0])

    assert repr(S0) == "Sum([FloatValue(%s)])" % f_0
    assert S0.t == CONST
    assert repr(S1) == "Sum([Symbol('x', BASIS)])"
    assert repr(S4) == "Sum([FloatValue(-%s), Symbol('x', BASIS)])" % f_2
    assert repr(S9) == "Sum([Symbol('y', GEO)])"

    assert str(S2) == "(x + y)"
    assert str(S3) == "(x + x)"
    assert str(S5) == "(x + x-%s)" % f_2
    assert str(S6) == "(%s + x + x)" % f_1
    assert str(S7) == "(x + y-%s)" % f_2
    assert str(S8) == "(%s + x + x)" % f_3
    assert str(S9) == "y"

    assert S2 == S2
    assert S2 != S3
    assert S5 != S6
    assert S2 == S10

    assert S0.ops() == 0
    assert S1.ops() == 0
    assert S2.ops() == 1
    assert S3.ops() == 1
    assert S4.ops() == 1
    assert S5.ops() == 2
    assert S6.ops() == 2
    assert S7.ops() == 2
    assert S8.ops() == 2
    assert S9.ops() == 0

    # Test hash
    ll = [S2]
    d = {S2: 0}

    assert S2 in ll
    assert S2 in d
    assert S10 in ll
    assert S10 in d
