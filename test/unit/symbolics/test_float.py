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
set_float_formatting(FFC_PARAMETERS["precision"])


def testFloat():
    "Test simple FloatValue instance."
    f0 = FloatValue(1.5)
    f1 = FloatValue(-5)
    f2 = FloatValue(-1e-14)
    f3 = FloatValue(-1e-11)
    f4 = FloatValue(1.5)

    assert repr(f0) == "FloatValue(%s)" % format["float"](1.5)
    assert repr(f1) == "FloatValue(%s)" % format["float"](-5)
    # This depends on the chosen precision and certain rounding behaviour in code generation:
    #assert repr(f2) == "FloatValue(%s)" % format["float"](0.0)
    assert repr(f3) == "FloatValue(%s)" % format["float"](-1e-11)

    #assert f2.val == 0  # This test documents incorrect float behaviour of quadrature repsentation
    assert not f3.val == 0

    assert f0.ops() == 0
    assert f1.ops() == 0
    assert f2.ops() == 0
    assert f3.ops() == 0

    assert f0 == f4
    assert f1 != f3
    assert not f0 < f1
    #assert f2 > f3
    # ^^^ NB! The > operator of FloatValue is not a comparison of numbers,
    # it does something else entirely, and is affected by precision in indirect ways.

    # Test hash
    ll = [f0]
    d = {f0: 0}
    assert f0 in ll
    assert f0 in d
    assert f4 in ll
    assert f4 in d
    assert f1 not in ll
    assert f1 not in d
