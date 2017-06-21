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


def testReduceOperations():

    f_1 = format["float"](1)
    f_2 = format["float"](2)

    # Aux. variables
    f2 = FloatValue(2)
    f0_5 = FloatValue(0.5)
    f1 = FloatValue(1.0)
    fm1 = FloatValue(-1.0)

    x = Symbol("x", GEO)
    y = Symbol("y", GEO)
    z = Symbol("z", GEO)
    a = Symbol("a", GEO)
    b = Symbol("b", GEO)
    c = Symbol("c", GEO)
    d = Symbol("d", GEO)

    # Simple expand and reduce simple float and symbol objects
    fx2 = f2.expand()
    xx = x.expand()

    fr2 = fx2.reduce_ops()
    xr = xx.reduce_ops()

    assert f2 == fr2
    assert x == xr

    # Test product
    p0 = f2 * x
    p1 = y * x
    p2 = x * f2 / y
    p3 = x * Sum([x, y])

    px0 = p0.expand()
    px1 = p1.expand()

    pr0 = px0.reduce_ops()
    pr1 = px1.reduce_ops()

    assert p0 == pr0
    assert p1 == pr1

    # Test fraction
    F0 = Fraction(p0, y)
    F1 = Fraction(x, p0)
    F2 = Fraction(p0, p1)
    F3 = Fraction(Sum([x * x, x * y]), y)
    F4 = Fraction(Sum([f2 * x, x * y]), a)

    Fx0 = F0.expand()
    Fx1 = F1.expand()
    Fx2 = F2.expand()
    Fx3 = F3.expand()
    Fx4 = F4.expand()

    Fr0 = Fx0.reduce_ops()
    Fr1 = Fx1.reduce_ops()
    Fr2 = Fx2.reduce_ops()
    Fr3 = Fx3.reduce_ops()
    Fr4 = Fx4.reduce_ops()

    assert Fr0 == F0
    assert Fr1 == f0_5
    assert Fr2 == Fraction(f2, y)
    assert str(Fr3) == "x*(%s + x/y)" % f_1
    assert str(Fr4) == "x*(%s + y)/a" % f_2

    # Test sum
    # TODO: Here we might have to add additional tests
    S0 = Sum([x, y])
    S1 = Sum([p0, p1])
    S2 = Sum([x, p1])
    S3 = Sum([p0, f2 * y])
    S4 = Sum([f2 * p1, z * p1])
    S5 = Sum([x, x * x, x * x * x])
    S6 = Sum([a * x * x, b * x * x * x, c * x * x, d * x * x * x])
    S7 = Sum([p0, p1, x * x, f2 * z, y * z])
    S8 = Sum([a * y, b * y, x * x * x * y, x * x * x * z])
    S9 = Sum([a * y, b * y, c * y, x * x * x * y, f2 * x * x, x * x * x * z])
    S10 = Sum([f2 * x * x * y, x * x * y * z])
    S11 = Sum([f2 * x * x * y * y, x * x * y * y * z])
    S12 = Sum([f2 * x * x * y * y, x * x * y * y * z, a * z, b * z, c * z])
    S13 = Sum([Fraction(f1, x), Fraction(f1, y)])
    S14 = Sum([Fraction(fm1, x), Fraction(fm1, y)])
    S15 = Sum([Fraction(f2, x), Fraction(f2, x)])
    S16 = Sum([Fraction(f2 * x, y * z), Fraction(f0_5, y * z)])
    S17 = Sum([(f2 * x * y) / a, (x * y * z) / b])
    S18 = Sum([(x * y) / a, (x * z) / a, f2 / a, (f2 * x * y) / a])
    S19 = Sum([(f2 * x) / a, (x * y) / a, z * x])
    S20 = Product([Sum([x, y]), Fraction(a, b), Fraction(Product([c, d]), z)])
    S21 = Sum([a * x, b * x, c * x, x * y, x * z, f2 * y, a * y, b * y, f2 * z, a * z, b * z])
    S22 = Sum([FloatValue(0.5) * x / y, FloatValue(-0.5) * x / y])
    S23 = Sum([x * y * z, x * y * y * y * z * z * z, y * y * y * z * z * z * z, z * z * z * z * z])

    Sx0 = S0.expand()
    Sx1 = S1.expand()
    Sx2 = S2.expand()
    Sx3 = S3.expand()
    Sx4 = S4.expand()
    Sx5 = S5.expand()
    Sx6 = S6.expand()
    Sx7 = S7.expand()
    Sx8 = S8.expand()
    Sx9 = S9.expand()
    Sx10 = S10.expand()
    Sx11 = S11.expand()
    Sx12 = S12.expand()
    Sx13 = S13.expand()
    Sx14 = S14.expand()
    Sx15 = S15.expand()
    Sx16 = S16.expand()
    Sx17 = S17.expand()
    Sx18 = S18.expand()
    Sx19 = S19.expand()
    Sx20 = S20.expand()
    Sx21 = S21.expand()
    Sx22 = S22.expand()
    Sx23 = S23.expand()

    Sr0 = Sx0.reduce_ops()
    Sr1 = Sx1.reduce_ops()
    Sr2 = Sx2.reduce_ops()
    Sr3 = Sx3.reduce_ops()
    Sr4 = Sx4.reduce_ops()
    Sr5 = Sx5.reduce_ops()
    Sr6 = Sx6.reduce_ops()
    Sr7 = Sx7.reduce_ops()
    Sr8 = Sx8.reduce_ops()
    Sr9 = Sx9.reduce_ops()
    Sr10 = Sx10.reduce_ops()
    Sr11 = Sx11.reduce_ops()
    Sr12 = Sx12.reduce_ops()
    Sr13 = Sx13.reduce_ops()
    Sr14 = Sx14.reduce_ops()
    Sr15 = Sx15.reduce_ops()
    Sr16 = Sx16.reduce_ops()
    Sr17 = Sx17.reduce_ops()
    Sr18 = Sx18.reduce_ops()
    Sr19 = Sx19.reduce_ops()
    Sr20 = Sx20.reduce_ops()
    Sr21 = Sx21.reduce_ops()
    Sr22 = Sx22.reduce_ops()
    Sr23 = Sx23.reduce_ops()

    assert Sr0 == S0
    assert str(Sr1) == "x*(%s + y)" % f_2
    # TODO: Should this be (x + x*y)?
    assert str(Sr2) == "x*(%s + y)" % f_1
    # assert str(Sr2), "(x + x*y)")
    assert str(Sr3) == "%s*(x + y)" % f_2
    assert str(Sr4) == "x*y*(%s + z)" % f_2
    assert str(Sr5) == "x*(%s + x*(%s + x))" % (f_1, f_1)
    assert str(Sr6) == "x*x*(a + c + x*(b + d))"
    assert str(Sr7) == "(x*(%s + x + y) + z*(%s + y))" % (f_2, f_2)
    assert str(Sr8) == "(x*x*x*(y + z) + y*(a + b))"
    assert str(Sr9) == "(x*x*(%s + x*(y + z)) + y*(a + b + c))" % f_2
    assert str(Sr10) == "x*x*y*(%s + z)" % f_2
    assert str(Sr11) == "x*x*y*y*(%s + z)" % f_2
    assert str(Sr12) == "(x*x*y*y*(%s + z) + z*(a + b + c))" % f_2
    assert str(Sr13) == "(%s/x + %s/y)" % (f_1, f_1)
    assert str(Sr14) == "(-%s/x-%s/y)" % (f_1, f_1)
    assert str(Sr15) == "%s/x" % format["float"](4)
    assert str(Sr16) == "(%s + %s*x)/(y*z)" % (format["float"](0.5), f_2)
    assert str(Sr17) == "x*y*(%s/a + z/b)" % f_2
    assert str(Sr18) == "(%s + x*(z + %s*y))/a" % (f_2, format["float"](3))
    assert str(Sr19) == "x*(z + (%s + y)/a)" % f_2
    assert str(Sr20) == "a*c*d*(x + y)/(b*z)"
    assert str(Sr21) == "(x*(a + b + c + y + z) + y*(%s + a + b) + z*(%s + a + b))" % (f_2, f_2)
    assert str(Sr22) == "%s" % format["float"](0)
    assert str(Sr23) == "(x*y*z + z*z*z*(y*y*y*(x + z) + z*z))"

    assert S0.ops() == 1
    assert Sr0.ops() == 1
    assert S1.ops() == 3
    assert Sr1.ops() == 2
    assert S2.ops() == 2
    assert Sr2.ops() == 2
    assert S3.ops() == 3
    assert Sr3.ops() == 2
    assert S4.ops() == 5
    assert Sr4.ops() == 3
    assert S5.ops() == 5
    assert Sr5.ops() == 4
    assert S6.ops() == 13
    assert Sr6.ops() == 6
    assert S7.ops() == 9
    assert Sr7.ops() == 6
    assert S8.ops() == 11
    assert Sr8.ops() == 7
    assert S9.ops() == 16
    assert Sr9.ops() == 9
    assert S10.ops() == 7
    assert Sr10.ops() == 4
    assert S11.ops() == 9
    assert Sr11.ops() == 5
    assert S12.ops() == 15
    assert Sr12.ops() == 9
    assert S13.ops() == 3
    assert Sr13.ops() == 3
    assert S14.ops() == 3
    assert Sr14.ops() == 3
    assert S15.ops() == 3
    assert Sr15.ops() == 1
    assert S16.ops() == 6
    assert Sr16.ops() == 4
    assert S17.ops() == 7
    assert Sr17.ops() == 5
    assert S18.ops() == 11
    assert Sr18.ops() == 5
    assert S19.ops() == 7
    assert Sr19.ops() == 4
    assert S20.ops() == 6
    assert Sr20.ops() == 6
    assert S21.ops() == 21
    assert Sr21.ops() == 13
    assert S23.ops() == 21
    assert Sr23.ops() == 12
