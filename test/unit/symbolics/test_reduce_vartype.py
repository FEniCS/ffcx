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


def testReduceVarType():
    f1 = FloatValue(1)
    f2 = FloatValue(2)
    f3 = FloatValue(3)
    f5 = FloatValue(5)
    fm4 = FloatValue(-4)

    B0 = Symbol("B0", BASIS)
    B1 = Symbol("B1", BASIS)
    Bm4 = Product([fm4, B1])
    B5 = Product([f5, B0])

    I0 = Symbol("I0", IP)
    I1 = Symbol("I1", IP)
    I2 = Symbol("I2", IP)
    I5 = Product([f5, I0])

    G0 = Symbol("G0", GEO)
    G1 = Symbol("G1", GEO)
    G2 = Symbol("G2", GEO)
    G3 = Product([f3, G0])

    C0 = Symbol("C0", CONST)
    C2 = Product([f2, C0])

    p0 = Product([B0, I5])
    p1 = Product([B0, B1])

    S0 = Sum([B0, I5])
    S1 = Sum([p0, p1])
    S2 = Sum([B0, B1])
    S3 = Sum([B0, p0])
    S4 = Sum([f5, p0])
    S5 = Sum([I0, G0])

    F0 = Fraction(B0, I5).expand()
    F1 = Fraction(p1, I5).expand()
    F2 = Fraction(G3, S2).expand()
    F3 = Fraction(G3, S3).expand()
    F4 = Fraction(I1, Sum([I1, I0]))
    F5 = Fraction(S5, I1)
    F6 = Fraction(I0,
                  Sum([
                      Fraction(Sum([I0, I1]), Sum([G0, G1])),
                      Fraction(Sum([I1, I2]), Sum([G1, G2])),
                  ]))

    r0 = B0.reduce_vartype(BASIS)
    r1 = B0.reduce_vartype(CONST)

    rp0 = p0.reduce_vartype(BASIS)
    rp1 = p0.reduce_vartype(IP)
    rp2 = p1.reduce_vartype(BASIS)
    rp3 = p1.reduce_vartype(GEO)

    rs0 = S0.reduce_vartype(BASIS)
    rs1 = S0.reduce_vartype(IP)
    rs2 = S1.reduce_vartype(BASIS)
    rs3 = S4.reduce_vartype(BASIS)
    rs4 = S4.reduce_vartype(CONST)

    rf0 = F0.reduce_vartype(BASIS)
    rf1 = F1.reduce_vartype(BASIS)
    rf2 = F0.reduce_vartype(IP)
    rf3 = F2.reduce_vartype(BASIS)
    rf4 = F3.reduce_vartype(BASIS)
    rf5 = F4.reduce_vartype(IP)
    rf6 = F5.reduce_vartype(IP)
    rf7 = F6.reduce_vartype(IP)

    assert [(B0, f1)] == r0
    assert [((), B0)] == r1

    assert [(B0, I5)] == rp0
    assert [(I0, B5)] == rp1
    assert [(p1, f1)] == rp2
    assert [((), p1)] == rp3

    assert ((), I5) == rs0[0]
    assert (B0, f1) == rs0[1]
    assert (I0, f5) == rs1[1]
    assert ((), B0) == rs1[0]
    assert (Product([B0, B1]), f1) == rs2[1]
    assert (B0, I5) == rs2[0]
    assert ((), f5) == rs3[0]
    assert (B0, I5) == rs3[1]
    assert (f5, Sum([f1, Product([B0, I0])])) == rs4[0]

    assert [(B0, Fraction(FloatValue(0.2), I0))] == rf0
    assert [(Product([B0, B1]), Fraction(FloatValue(0.2), I0))] == rf1
    assert [(Fraction(f1, I0), Product([FloatValue(0.2), B0]))] == rf2
    assert [(Fraction(f1, S2), G3)] == rf3
    assert [(Fraction(f1, B0), Fraction(G3, Sum([I5, f1])))] == rf4
    assert F4 == rf5[0][0]
    assert FloatValue(1) == rf5[0][1]
    assert Fraction(I0, I1) == rf6[1][0]
    assert f1 == rf6[1][1]
    assert Fraction(f1, I1) == rf6[0][0]
    assert G0 == rf6[0][1]
    assert F6 == rf7[0][0]
    assert f1 == rf7[0][1]

    expr = Sum([Symbol('W1', GEO), Fraction(Symbol('det', GEO), Sum([Symbol('F0', IP), Symbol('K_11', GEO)]))])
    red = expr.expand().reduce_vartype(IP)
    vals = []
    for ip in red:
        ip_dec, geo = ip
        if ip_dec and geo:
            vals.append(Product([ip_dec, geo]))
        elif geo:
            vals.append(geo)
        elif ip_dec:
            vals.append(ip_dec)
    comb = Sum(vals).expand()
    K_11 = 1.4
    F0 = 1.5
    W1 = 1.9
    det = 2.1
    assert round(eval(str(expr)) - eval(str(comb)), 10) == 0.0
