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


def testSymbol():
    "Test simple symbol instance."

    s0 = Symbol("x", BASIS)
    s1 = Symbol("y", IP)
    s2 = Symbol("z", GEO)
    s3 = Symbol("z", GEO)
    s4 = Symbol("z", IP)

    assert repr(s0) == "Symbol('x', BASIS)"
    assert repr(s1) == "Symbol('y', IP)"
    assert repr(s2) == "Symbol('z', GEO)"
    assert repr(s4) == "Symbol('z', IP)"

    assert s2 == s3
    assert (s2 == s1) is False
    assert (s2 == s4) is False
    assert (s2 != s3) is False
    assert s2 != s1

    assert s0 < s1
    assert s4 > s1

    assert s0.ops() == 0
    assert s1.ops() == 0
    assert s2.ops() == 0
    assert s3.ops() == 0
    assert s4.ops() == 0

    # Test hash
    ll = [s0]
    d = {s0: 0}
    s5 = Symbol('x', BASIS)

    assert s0 in ll
    assert s0 in d
    assert s5 in ll
    assert s5 in d
