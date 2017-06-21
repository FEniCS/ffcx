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


def test_Elasticity_Term():
    # expr:  0.25*W1*det*(FE0_C2_D001[0][j]*FE0_C2_D001[0][k]*Jinv_00*Jinv_21 + FE0_C2_D001[0][j]*FE0_C2_D001[0][k]*Jinv_00*Jinv_21)
    expr = Product([
        FloatValue(0.25), Symbol('W1', GEO), Symbol('det', GEO),
        Sum([Product([Symbol('FE0_C2_D001_0_j', BASIS), Symbol('FE0_C2_D001_0_k', BASIS),
                      Symbol('Jinv_00', GEO), Symbol('Jinv_21', GEO)]),
             Product([Symbol('FE0_C2_D001_0_j', BASIS), Symbol('FE0_C2_D001_0_k', BASIS),
                      Symbol('Jinv_00', GEO), Symbol('Jinv_21', GEO)])
             ])
    ])

    expr_exp = expr.expand()
    expr_red = expr_exp.reduce_ops()

    det, W1, Jinv_00, Jinv_21, FE0_C2_D001_0_j, FE0_C2_D001_0_k = [0.123 + i for i in range(6)]
    assert round(eval(str(expr)) - eval(str(expr_exp)), 10) == 0.0
    assert round(eval(str(expr)) - eval(str(expr_red)), 10) == 0.0
    assert expr.ops() == 10
    assert expr_exp.ops() == 6
    assert expr_red.ops() == 6

    # Generate code
    ip_consts = {}
    geo_consts = {}
    trans_set = set()

    start = time.time()
    opt_code = optimise_code(expr, ip_consts, geo_consts, trans_set)

    G = [eval(str(list(geo_consts.items())[0][0]))]
    assert round(eval(str(expr)) - eval(str(opt_code))) == 0.0
