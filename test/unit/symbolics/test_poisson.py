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
from ffc.quadrature.reduce_operations import operation_count, expand_operations, reduce_operations
from ffc.quadrature.symbolics import *
from ffc.quadrature.cpp import format, set_float_formatting
from ffc.parameters import FFC_PARAMETERS
set_float_formatting(FFC_PARAMETERS['precision'])


def testPoisson():
    poisson = """((Jinv_00*FE0_D10_ip_j + Jinv_10*FE0_D01_ip_j)*(Jinv_00*FE0_D10_ip_k + Jinv_10*FE0_D01_ip_k) + (Jinv_01*FE0_D10_ip_j + Jinv_11*FE0_D01_ip_j)*(Jinv_01*FE0_D10_ip_k + Jinv_11*FE0_D01_ip_k))*W4_ip*det"""

    expr = Product([
        Sum([
            Product([
                Sum([
                    Product([Symbol("Jinv_00", GEO), Symbol("FE0_D10_ip_j", BASIS)]),
                    Product([Symbol("Jinv_10", GEO), Symbol("FE0_D01_ip_j", BASIS)])
                ]),
                Sum([
                    Product([Symbol("Jinv_00", GEO), Symbol("FE0_D10_ip_k", BASIS)]),
                    Product([Symbol("Jinv_10", GEO), Symbol("FE0_D01_ip_k", BASIS)])
                ])
            ]),
            Product([
                Sum([
                    Product([Symbol("Jinv_01", GEO), Symbol("FE0_D10_ip_j", BASIS)]),
                    Product([Symbol("Jinv_11", GEO), Symbol("FE0_D01_ip_j", BASIS)])
                ]),
                Sum([
                    Product([Symbol("Jinv_01", GEO), Symbol("FE0_D10_ip_k", BASIS)]),
                    Product([Symbol("Jinv_11", GEO), Symbol("FE0_D01_ip_k", BASIS)])
                ])
            ])
        ]),
        Symbol("W4_ip", IP),
        Symbol("det", GEO)
    ])

    expr_exp = expr.expand()
    poisson_exp = expand_operations(poisson, format)

    expr_red = expr_exp.reduce_ops()
    poisson_red = reduce_operations(poisson, format)

    poisson_exp_ops = operation_count(poisson_exp, format)
    poisson_red_ops = operation_count(poisson_red, format)

    Jinv_00, Jinv_01, Jinv_10, Jinv_11, W4_ip, det = (1.1, 1.5, -4.3, 1.7, 11, 52.3)
    FE0_D01_ip_j, FE0_D10_ip_j, FE0_D01_ip_k, FE0_D10_ip_k = (1.12, 5.7, -9.3, 7.4)
    assert round(eval(str(expr)) - eval(str(expr_exp)), 10) == 0.0
    assert round(eval(str(expr)) - eval(str(expr_red)), 10) == 0.0
    assert round(eval(str(expr)) - eval(str(poisson)), 10) == 0.0
    assert round(eval(str(expr)) - eval(str(poisson_exp)), 10) == 0.0
    assert round(eval(str(expr)) - eval(str(poisson_red)), 10) == 0.0
    assert expr.ops() == 17
    assert poisson_exp_ops == 47
    assert expr_exp.ops() == 47
    assert poisson_red_ops == 23
    assert expr_red.ops() == 23
