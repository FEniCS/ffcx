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


def testDGElastoDyn():
    expr = Product([
        Sum([
            Symbol("F0", IP),
            Symbol("F1", IP)
        ]),
        Fraction(
            Symbol("w4", GEO),
            Symbol("w3", GEO)
        ),
        Fraction(
            Product([
                Symbol("w2", GEO),
                Symbol("w5", GEO)
            ]),
            Symbol("w6", GEO)
        )
    ])

    expr_exp = expr.expand()

    expr_red = expr_exp.reduce_ops()

    F0, F1, w2, w3, w4, w5, w6 = (3.12, -8.1, -45.3, 17.5, 2.2, 5.3, 9.145)
    assert round(eval(str(expr)) - eval(str(expr_exp)), 10) == 0
    assert round(eval(str(expr)) - eval(str(expr_red)), 10) == 0
    assert expr.ops() == 6
    assert expr_exp.ops() == 11
    assert expr_red.ops() == 6
