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
from ffc.log import push_level, pop_level, CRITICAL


def testFraction():
    "Test simple fraction instance."

    f0 = FloatValue(-2.0)
    f1 = FloatValue(3.0)
    f2 = FloatValue(0)
    s0 = Symbol("x", BASIS)
    s1 = Symbol("y", GEO)

    F0 = Fraction(f1, f0)
    F1 = Fraction(f2, f0)
    F2 = Fraction(s0, s1)
    F3 = Fraction(s0, f1)
    F4 = Fraction(f0, s1)
    F5 = Fraction(f2, s1)
    F6 = Fraction(s0, s1)

    # Silence output
    push_level(CRITICAL)
    with pytest.raises(Exception):
        Fraction(f0, f2)
    with pytest.raises(Exception):
        Fraction(s0, f2)
    pop_level()

    assert repr(F0) == "Fraction(FloatValue(%s), FloatValue(%s))"\
        % (format["float"](-1.5), format["float"](1))
    assert repr(F2) == "Fraction(Symbol('x', BASIS), Symbol('y', GEO))"

    assert str(F0) == "%s" % format["float"](-1.5)
    assert str(F1) == "%s" % format["float"](0)
    assert str(F2) == "x/y"
    assert str(F3) == "x/%s" % format["float"](3)
    assert str(F4) == "-%s/y" % format["float"](2)
    assert str(F5) == "%s" % format["float"](0)

    assert F2 == F2
    assert F2 != F3
    assert F5 != F4
    assert F2 == F6

    assert F0.ops() == 0
    assert F1.ops() == 0
    assert F2.ops() == 1
    assert F3.ops() == 1
    assert F4.ops() == 1
    assert F5.ops() == 0

    # Test hash
    ll = [F2]
    d = {F2: 0}

    assert F2 in ll
    assert F2 in d
    assert F6 in ll
    assert F6 in d
