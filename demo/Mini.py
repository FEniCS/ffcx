# Copyright (C) 2010 Marie E. Rognes
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
# Illustration of vector sum of elements (EnrichedElement): The
# bilinear form a(u, v) for the Stokes equations using a mixed
# formulation involving the Mini element. The velocity element is
# composed of a P1 element augmented by the cubic bubble function.
from basix.ufl import BlockedElement, MixedElement, element, enriched_element
from ufl import TestFunctions, TrialFunctions, div, dx, grad, inner

P1 = element("Lagrange", "triangle", 1)
B = element("Bubble", "triangle", 3)
V = BlockedElement(enriched_element([P1, B]), shape=(2, ))
Q = element("CG", "triangle", 1)

Mini = MixedElement([V, Q])
(u, p) = TrialFunctions(Mini)
(v, q) = TestFunctions(Mini)

a = (inner(grad(u), grad(v)) - div(v) * p + div(u) * q) * dx
