# Copyright (C) 2020 Matthew Scroggs
#
# This file is part of FFCX.
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
#
"Unit tests for FFCX"


import numpy
import pytest

from ffcx.fiatinterface import create_element
from ffcx.ir.representation import _get_coeffs_from_fiat, _get_dmats_from_fiat
from ufl import FiniteElement


def test_rt_triangle():
    P = create_element(FiniteElement("RT", "triangle", 1))
    co = _get_coeffs_from_fiat(P)

    points = []
    for i in P.dual_basis():
        points += list(i.pt_dict.keys())

    tab = P.tabulate(0, points)[0,0]


    from IPython import embed; embed()


def test_rtce_quad():
    P = create_element(FiniteElement("RTCE", "quadrilateral", 1))
    from IPython import embed; embed()

    for e in P.elements():
        _get_dmats_from_fiat(e)

    co = _get_coeffs_from_fiat(P)

    points = []
    for i in P.dual_basis():
        points += list(i.pt_dict.keys())

    tab = P.tabulate(0, points)[0,0]


    from IPython import embed; embed()
