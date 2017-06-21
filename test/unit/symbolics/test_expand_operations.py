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


def testExpandOperations():
    f0 = FloatValue(-1)
    f1 = FloatValue(2)
    f2 = FloatValue(1)
    sx = Symbol("x", GEO)
    sy = Symbol("y", GEO)
    sz = Symbol("z", GEO)
    s0 = Product([FloatValue(-1), Symbol("x", GEO)])
    s1 = Symbol("y", GEO)
    s2 = Product([FloatValue(5), Symbol("z", IP)])
    s3 = Product([FloatValue(-4), Symbol("z", GEO)])

    # Random variable values
    x = 2.2
    y = -0.2
    z = 1.1

    # Aux. expressions
    P0 = Product([s2, s1])
    P1 = Product([P0, s0])
    P2 = Product([P1, s1, P0])
    P3 = Product([P1, P2])

    S0 = Sum([s2, s1])
    S1 = Sum([S0, s0])
    S2 = Sum([S1, s1, S0])
    S3 = Sum([S1, S2])

    F0 = Fraction(s2, s1)
    F1 = Fraction(F0, s0)
    F2 = Fraction(F1, F0)
    F3 = Fraction(F1, F2)

    # Special fractions
    F4 = Fraction(P0, F0)
    F5 = Fraction(Fraction(s0, P0), P0)
    F6 = Fraction(Fraction(Fraction(s1, s0), Fraction(s1, s2)), Fraction(Fraction(s2, s0), Fraction(s1, s0)))
    F7 = Fraction(s1, Product([s1, Symbol("x", GEO)]))
    F8 = Fraction(Sum([sx, Fraction(sy, sx)]), FloatValue(2))

    F4x = F4.expand()
    F5x = F5.expand()
    F6x = F6.expand()
    F7x = F7.expand()
    F8x = F8.expand()

    assert round(eval(str(F4)) - eval(str(F4x)), 10) == 0.0
    assert round(eval(str(F5)) - eval(str(F5x)), 10) == 0.0
    assert round(eval(str(F6)) - eval(str(F6x)), 10) == 0.0
    assert round(eval(str(F7)) - eval(str(F7x)), 10) == 0.0
    assert round(eval(str(F8)) - eval(str(F8x)), 10) == 0.0

    assert F4.ops() == 5
    assert F4x.ops() == 1
    assert F5.ops() == 6
    assert F5x.ops() == 5
    assert F6.ops() == 9
    assert F6x.ops() == 1
    assert F7.ops() == 2
    assert F7x.ops() == 1
    assert F8.ops() == 3
    assert F8x.ops() == 4

    # Expressions that should be expanded
    e0 = Product([P3, F2])
    e1 = Product([S3, P2])
    e2 = Product([F3, S1])

    e3 = Sum([P3, F2])
    e4 = Sum([S3, P2])
    e5 = Sum([F3, S1])

    e6 = Fraction(P3, F2)
    e7 = Fraction(S3, P2)
    e8 = Fraction(F3, S1)
    e9 = Fraction(S0, s0)

    e0x = e0.expand()
    e1x = e1.expand()
    e2x = e2.expand()
    e3x = e3.expand()
    e4x = e4.expand()
    e5x = e5.expand()
    e6x = e6.expand()
    e7x = e7.expand()
    e8x = e8.expand()
    e9x = e9.expand()

    assert round(eval(str(e0)) - eval(str(e0x)), 10) == 0.0
    assert round(eval(str(e1)) - eval(str(e1x)), 10) == 0.0
    assert round(eval(str(e2)) - eval(str(e2x)), 10) == 0.0
    assert round(eval(str(e3)) - eval(str(e3x)), 10) == 0.0
    assert round(eval(str(e4)) - eval(str(e4x)), 10) == 0.0
    assert round(eval(str(e5)) - eval(str(e5x)), 10) == 0.0
    assert round(eval(str(e6)) - eval(str(e6x)), 10) == 0.0
    assert round(eval(str(e7)) - eval(str(e7x)), 10) == 0.0
    assert round(eval(str(e8)) - eval(str(e8x)), 10) == 0.0
    assert round(eval(str(e9)) - eval(str(e9x)), 10) == 0.0

    assert e0.ops() == 16
    assert e0x.ops() == 8
    assert e1.ops() == 18
    assert e1x.ops() == 23
    assert e2.ops() == 14
    assert e2x.ops() == 9
    assert e3.ops() == 16
    assert e3x.ops() == 11
    assert e4.ops() == 18
    assert e4x.ops() == 12
    assert e5.ops() == 14
    assert e5x.ops() == 6
    assert e6.ops() == 16
    assert e6x.ops() == 10
    assert e7.ops() == 18
    assert e7x.ops() == 17
    assert e8.ops() == 14
    assert e8x.ops() == 8
    assert e9.ops() == 3
    assert e9x.ops() == 4

    # More expressions (from old expand tests)
    PF = Product([F0, F1])
    E0 = Product([s1, f0, S1])
    E1 = Sum([P0, E0])
    E2 = Fraction(Sum([Product([f1])]), f2)
    E3 = Sum([F0, F0])
    E4 = Product([Sum([Product([sx, Sum([sy, Product([Sum([sy, Product([sy, sz]), sy])]), sy])]),
                       Product([sx, Sum([Product([sy, sz]), sy])])])])
    P4 = Product([s1, Sum([s0, s1])])
    P5 = Product([s0, E0])
    P6 = Product([s1])
    S4 = Sum([s1])

    # Create 'real' term that caused me trouble
    P00 = Product([Symbol("Jinv_00", GEO)] * 2)
    P01 = Product([Symbol("Jinv_01", GEO)] * 2)
    P20 = Product([Symbol("Jinv_00", GEO),
                   Product([f1, Symbol("Jinv_20", GEO)])])
    P21 = Product([Symbol("Jinv_01", GEO),
                   Product([f1, Symbol("Jinv_21", GEO)])])
    PS0 = Product([Symbol("Jinv_22", GEO),
                   Sum([P00, P01])])
    PS1 = Product([Product([f0, Symbol("Jinv_02", GEO)]),
                   Sum([P20, P21])])
    SP = Sum([PS0, PS1])

    PFx = PF.expand()
    E0x = E0.expand()
    E1x = E1.expand()
    E2x = E2.expand()
    E3x = E3.expand()
    E4x = E4.expand()
    P4x = P4.expand()
    P5x = P5.expand()
    P6x = P6.expand()
    S4x = S4.expand()
    SPx = SP.expand()

    Jinv_00, Jinv_01, Jinv_10, Jinv_02, Jinv_20, Jinv_22, Jinv_21, W1, det = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    assert round(eval(str(SP)) - eval(str(SPx)), 10) == 0.0
    assert round(eval(str(E0)) - eval(str(E0x)), 10) == 0.0
    assert round(eval(str(E1)) - eval(str(E1x)), 10) == 0.0
    assert round(eval(str(E2)) - eval(str(E2x)), 10) == 0.0
    assert round(eval(str(E3)) - eval(str(E3x)), 10) == 0.0
    assert round(eval(str(E4)) - eval(str(E4x)), 10) == 0.0
    assert round(eval(str(SP)) - eval(str(SPx)), 10) == 0.0
    assert round(eval(str(P4)) - eval(str(P4x)), 10) == 0.0
    assert round(eval(str(P5)) - eval(str(P5x)), 10) == 0.0
    assert P6x == s1
    assert S4x == s1
    assert PF.ops() == 6
    assert PFx.ops() == 5
    assert E0.ops() == 4
    assert E0x.ops() == 6
    assert E1.ops() == 7
    assert E1x.ops() == 3
    assert E2.ops() == 1
    assert E2x.ops() == 0
    assert E3.ops() == 5
    assert E3x.ops() == 5
    assert E4.ops() == 10
    assert E4x.ops() == 6
    assert SP.ops() == 11
    assert SPx.ops() == 13
    assert P4.ops() == 2
    assert P4x.ops() == 3
    assert P5.ops() == 5
    assert P5x.ops() == 9
