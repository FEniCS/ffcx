# Copyright (C) 2015 Marie E. Rognes
#
# This file is part of FFCx.
#
# FFCx is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFCx is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFCx. If not, see <http://www.gnu.org/licenses/>.
"""Trace element demo."""

import basix.ufl
from ufl import FunctionSpace, Mesh, TestFunction, avg, dS, ds

element = basix.ufl.element("HDiv Trace", "triangle", 0)
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, element)
v = TestFunction(space)
L = v * ds + avg(v) * dS
