# Copyright (C) 2004-2007 Anders Logg
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
# The bilinear form a(u, v) and linear form L(v) for
# Poisson's equation.
#
# Compile this form with FFCx: ffcx Poisson.ufl
from ufl import (Coefficient, Constant, Mesh, TestFunction,
                 TrialFunction, dx, grad, inner)
import basix.ufl

mesh = Mesh(basix.ufl.element('P', "triangle", 2, shape=(2, )))

e = basix.ufl.element("Lagrange", "triangle", 2)

u = TrialFunction(e)
v = TestFunction(e)
f = Coefficient(e)

kappa1 = Constant(mesh.ufl_cell(), shape=(2, 2))
kappa2 = Constant(mesh.ufl_cell(), shape=(2, 2))

a = inner(kappa1, kappa2) * inner(grad(u), grad(v)) * dx
L = f * v * dx
