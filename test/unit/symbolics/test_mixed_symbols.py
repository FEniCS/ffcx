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


def testMixedSymbols():

    f_0 = format["float"](0)
    f_2 = format["float"](2)
    f_3 = format["float"](3)
    f_4 = format["float"](4)
    f_6 = format["float"](6)

    f0 = FloatValue(-2.0)
    f1 = FloatValue(3.0)
    f2 = FloatValue(0)

    s0 = Symbol("x", BASIS)
    s1 = Symbol("y", GEO)
    s2 = Symbol("z", GEO)

    p0 = Product([s0, s1])
    p1 = Product([f1, s0, s1])
    p2 = Product([s0, f2, s2])
    p3 = Product([s0, f0, s1, f1, s2])

    S0 = Sum([s0, s1])
    S1 = Sum([s0, s0])
    S2 = Sum([f0, s0])
    S3 = Sum([s0, f0, s0])

    F0 = Fraction(f1, f0)
    F1 = Fraction(s0, s1)
    F2 = Fraction(s0, f1)
    F3 = Fraction(f0, s1)

    x = 1.2
    y = 2.36
    z = 6.75
    # Mixed products
    mpp0 = Product([p0, s0])
    mpp1 = Product([p1, p0])
    mpp2 = Product([p2, p3])
    mpp3 = Product([p1, mpp1])

    mps0 = Product([S0, s0])
    mps1 = Product([S1, S0])
    mps2 = Product([S2, S3])
    mps3 = Product([S1, mps1])

    mpf0 = Product([F1, s0])
    mpf1 = Product([F1, F2])
    mpf2 = Product([F2, F3])
    mpf3 = Product([F1, mpf1])

    assert round(eval(str(mpp0)) - eval(str(p0)) * eval(str(s0)), 10) == 0.0
    assert round(eval(str(mpp1)) - eval(str(p1)) * eval(str(p0)), 10) == 0.0
    assert round(eval(str(mpp2)) - eval(str(p2)) * eval(str(p3)), 10) == 0.0
    assert round(eval(str(mpp3)) - eval(str(p1)) * eval(str(mpp1)), 10) == 0.0

    assert round(eval(str(mps0)) - eval(str(S0)) * eval(str(s0)), 10) == 0.0
    assert round(eval(str(mps1)) - eval(str(S1)) * eval(str(S0)), 10) == 0.0
    assert round(eval(str(mps2)) - eval(str(S2)) * eval(str(S3)), 10) == 0.0
    assert round(eval(str(mps3)) - eval(str(S1)) * eval(str(mps1)), 10) == 0.0

    assert round(eval(str(mpf0)) - eval(str(F1)) * eval(str(s0)), 10) == 0.0
    assert round(eval(str(mpf1)) - eval(str(F1)) * eval(str(F2)), 10) == 0.0
    assert round(eval(str(mpf2)) - eval(str(F2)) * eval(str(F3)), 10) == 0.0
    assert round(eval(str(mpf3)) - eval(str(F1)) * eval(str(mpf1)), 10) == 0.0

    assert mpp0.ops() == 2
    assert mpp1.ops() == 4
    assert mpp2.ops() == 0
    assert mpp3.ops() == 6

    assert mps0.ops() == 2
    assert mps1.ops() == 3
    assert mps2.ops() == 4
    assert mps3.ops() == 5

    assert mpf0.ops() == 2
    assert mpf1.ops() == 3
    assert mpf2.ops() == 3
    assert mpf3.ops() == 5

    assert str(mpp0) == 'x*x*y'
    assert str(mpp1) == '%s*x*x*y*y' % f_3
    assert str(mpp2) == '%s' % f_0
    assert str(mpp3) == '%s*x*x*x*y*y*y' % format["float"](9)
    assert str(mps0) == 'x*(x + y)'
    assert str(mps1) == '(x + x)*(x + y)'
    assert str(mps2) == '(x + x-%s)*(x-%s)' % (f_2, f_2)
    assert str(mps3) == '(x + x)*(x + x)*(x + y)'
    assert str(mpf0) == 'x*x/y'
    assert str(mpf1) == 'x/%s*x/y' % f_3
    assert str(mpf2) == '-%s/y*x/%s' % (f_2, f_3)
    assert str(mpf3) == 'x/%s*x/y*x/y' % f_3

    # Mixed sums
    msp0 = Sum([p0, s0])
    msp1 = Sum([p1, p0])
    msp2 = Sum([p2, p3])
    msp3 = Sum([p1, msp1])
    msp4 = Sum([f2, f2])

    mss0 = Sum([S0, s0])
    mss1 = Sum([S1, S0])
    mss2 = Sum([S2, S3])
    mss3 = Sum([S1, mps1])

    msf0 = Sum([F1, s0])
    msf1 = Sum([F1, F2])
    msf2 = Sum([F2, F3])
    msf3 = Sum([F1, msf1])

    assert round(eval(str(msp0)) - (eval(str(p0)) + eval(str(s0))), 10) == 0.0
    assert round(eval(str(msp1)) - (eval(str(p1)) + eval(str(p0))), 10) == 0.0
    assert round(eval(str(msp2)) - (eval(str(p2)) + eval(str(p3))), 10) == 0.0
    assert round(eval(str(msp3)) - (eval(str(p1)) + eval(str(msp1))), 10) == 0.0
    assert str(msp4) == '%s' % f_0

    assert round(eval(str(mss0)) - (eval(str(S0)) + eval(str(s0))), 10) == 0.0
    assert round(eval(str(mss1)) - (eval(str(S1)) + eval(str(S0))), 10) == 0.0
    assert round(eval(str(mss2)) - (eval(str(S2)) + eval(str(S3))), 10) == 0.0
    assert round(eval(str(mss3)) - (eval(str(S1)) + eval(str(mps1))), 10) == 0.0

    assert round(eval(str(msf0)) - (eval(str(F1)) + eval(str(s0))), 10) == 0.0
    assert round(eval(str(msf1)) - (eval(str(F1)) + eval(str(F2))), 10) == 0.0
    assert round(eval(str(msf2)) - (eval(str(F2)) + eval(str(F3))), 10) == 0.0
    assert round(eval(str(msf3)) - (eval(str(F1)) + eval(str(msf1))), 10) == 0.0

    assert msp0.ops() == 2
    assert msp1.ops() == 4
    assert msp2.ops() == 3
    assert msp3.ops() == 7

    assert mss0.ops() == 2
    assert mss1.ops() == 3
    assert mss2.ops() == 3
    assert mss3.ops() == 5

    assert msf0.ops() == 2
    assert msf1.ops() == 3
    assert msf2.ops() == 3
    assert msf3.ops() == 5

    assert str(msp0) == '(x + x*y)'
    assert str(msp1) == '(%s*x*y + x*y)' % f_3
    assert str(msp2) == '-%s*x*y*z' % f_6
    assert str(msp3) == '(%s*x*y + %s*x*y + x*y)' % (f_3, f_3)
    assert str(mss0) == '(x + x + y)'
    assert str(mss1) == '(x + x + x + y)'
    assert str(mss2) == '(x + x + x-%s)' % f_4
    assert str(mss3) == '(x + x + (x + x)*(x + y))'
    assert str(msf0) == '(x + x/y)'
    assert str(msf1) == '(x/%s + x/y)' % f_3
    assert str(msf2) == '(x/%s-%s/y)' % (f_3, f_2)
    assert str(msf3) == '(x/%s + x/y + x/y)' % f_3

    # Mixed fractions
    mfp0 = Fraction(p0, s0)
    mfp1 = Fraction(p1, p0)
    mfp2 = Fraction(p2, p3)
    mfp3 = Fraction(p1, mfp1)

    mfs0 = Fraction(S0, s0)
    mfs1 = Fraction(S1, S0)
    mfs2 = Fraction(S2, S3)
    mfs3 = Fraction(S1, mfs1)

    mff0 = Fraction(F1, s0)
    mff1 = Fraction(F1, F2)
    mff2 = Fraction(F2, F3)
    mff3 = Fraction(F1, mff1)

    assert round(eval(str(mfp0)) - eval(str(p0)) / eval(str(s0)), 10) == 0.0
    assert round(eval(str(mfp1)) - eval(str(p1)) / eval(str(p0)), 10) == 0.0
    assert round(eval(str(mfp2)) - eval(str(p2)) / eval(str(p3)), 10) == 0.0
    assert round(eval(str(mfp3)) - eval(str(p1)) / eval(str(mfp1)), 10) == 0.0

    assert round(eval(str(mfs0)) - eval(str(S0)) / eval(str(s0)), 10) == 0.0
    assert round(eval(str(mfs1)) - eval(str(S1)) / eval(str(S0)), 10) == 0.0
    assert round(eval(str(mfs2)) - eval(str(S2)) / eval(str(S3)), 10) == 0.0
    assert round(eval(str(mfs3)) - eval(str(S1)) / eval(str(mfs1)), 10) == 0.0

    assert round(eval(str(mff0)) - eval(str(F1)) / eval(str(s0)), 10) == 0.0
    assert round(eval(str(mff1)) - eval(str(F1)) / eval(str(F2)), 10) == 0.0
    assert round(eval(str(mff2)) - eval(str(F2)) / eval(str(F3)), 10) == 0.0
    assert round(eval(str(mff3)) - eval(str(F1)) / eval(str(mff1)), 10) == 0.0

    assert mfp0.ops() == 2
    assert mfp1.ops() == 4
    assert mfp2.ops() == 0
    assert mfp3.ops() == 7

    assert mfs0.ops() == 2
    assert mfs1.ops() == 3
    assert mfs2.ops() == 4
    assert mfs3.ops() == 5

    assert mff0.ops() == 2
    assert mff1.ops() == 3
    assert mff2.ops() == 3
    assert mff3.ops() == 5

    assert str(mfp0) == 'x*y/x'
    assert str(mfp1) == '%s*x*y/(x*y)' % f_3
    assert str(mfp2) == '%s' % f_0
    assert str(mfp3) == '%s*x*y/(%s*x*y/(x*y))' % (f_3, f_3)
    assert str(mfs0) == '(x + y)/x'
    assert str(mfs1) == '(x + x)/(x + y)'
    assert str(mfs2) == '(x-%s)/(x + x-%s)' % (f_2, f_2)
    assert str(mfs3) == '(x + x)/((x + x)/(x + y))'
    assert str(mff0) == '(x/y)/x'
    assert str(mff1) == '(x/y)/(x/%s)' % f_3
    assert str(mff2) == '(x/%s)/(-%s/y)' % (f_3, f_2)
    assert str(mff3) == '(x/y)/((x/y)/(x/%s))' % f_3

    # Use p1 as a base expression for Symbol
    s3 = Symbol(format["cos"](str(p1)), CONST, p1, 1)
    assert str(s3) == 'std::cos(%s*x*y)' % f_3
    assert s3.ops() == 3
