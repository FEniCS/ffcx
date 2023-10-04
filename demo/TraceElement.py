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
import basix.ufl
from ufl import TestFunction, avg, ds, dS

element = basix.ufl.element("HDiv Trace", "triangle", 0)
v = TestFunction(element)
L = v * ds + avg(v) * dS
