# Copyright (C) 2010 Kristian B. Oelgaard
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
#
# Test all algebra operators on Coefficients.
from ufl import (Coefficient, FiniteElement, acos, asin, atan, bessel_J,
                 bessel_Y, cos, dx, erf, exp, ln, sin, sqrt, tan, triangle)

element = FiniteElement("Lagrange", triangle, 1)

c0 = Coefficient(element)
c1 = Coefficient(element)

s0 = 3 * c0 - c1
p0 = c0 * c1
f0 = c0 / c1

integrand = sqrt(c0) + sqrt(s0) + sqrt(p0) + sqrt(f0)\
    + exp(c0) + exp(s0) + exp(p0) + exp(f0)\
    + ln(c0) + ln(s0) + ln(p0) + ln(f0)\
    + cos(c0) + cos(s0) + cos(p0) + cos(f0)\
    + sin(c0) + sin(s0) + sin(p0) + sin(f0)\
    + tan(c0) + tan(s0) + tan(p0) + tan(f0)\
    + acos(c0) + acos(s0) + acos(p0) + acos(f0)\
    + asin(c0) + asin(s0) + asin(p0) + asin(f0)\
    + atan(c0) + atan(s0) + atan(p0) + atan(f0)\
    + erf(c0) + erf(s0) + erf(p0) + erf(f0)\
    + bessel_J(1, c0) + bessel_J(1, s0) + bessel_J(0, p0) + bessel_J(0, f0)\
    + bessel_Y(1, c0) + bessel_Y(1, s0) + bessel_Y(0, p0) + bessel_Y(0, f0)

a = integrand * dx
