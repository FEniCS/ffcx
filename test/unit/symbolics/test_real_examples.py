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


def testRealExamples():
    p = Product([Symbol('FE0_C1_D01[ip][k]', BASIS),
                 Sum([
                     Symbol('Jinv_10', GEO),
                      Symbol('w[4][0]', GEO)
                     ]),
        Sum([
            Symbol('Jinv_10', GEO),
                      Symbol('w[4][0]', GEO)
        ])
    ])

    br = p.reduce_vartype(BASIS)
    be = p.expand().reduce_vartype(BASIS)

    if len(be) == 1:
        if be[0][0] == br[0]:
            if be[0][1] != br[1].expand():
                print("\nbe: ", repr(be[0][1]))
                print("\nbr: ", repr(br[1].expand()))
                print("\nbe: ", be[0][1])
                print("\nbr: ", br[1].expand())
                RuntimeError("Should not be here")
