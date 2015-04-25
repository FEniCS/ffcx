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
#
# First added:  2014-01-28
# Last changed: 2014-01-28

from ufl import FiniteElement, MixedElement
from .test_common import xcomb
__all__ = ["single_elements", "mixed_elements"]

# Elements, supported by FFC and FIAT, and their supported shape and orders
single_elements = [ {"family": "Lagrange",\
                      "shapes": ["interval", "triangle", "tetrahedron"],\
                      "orders": [1, 2, 3, 4]},\
                    {"family": "Discontinuous Lagrange",\
                      "shapes": ["interval", "triangle", "tetrahedron"],\
                      "orders": [0, 1, 2, 3, 4]},\
                    {"family": "Crouzeix-Raviart",\
                      "shapes": ["triangle", "tetrahedron"],\
                      "orders": [1]},\
                    {"family": "Raviart-Thomas",\
                      "shapes": ["triangle", "tetrahedron"],\
                      "orders": [1, 2, 3]},\
                    {"family": "Discontinuous Raviart-Thomas",\
                      "shapes": ["triangle", "tetrahedron"],\
                      "orders": [1, 2, 3]},\
                    {"family": "Brezzi-Douglas-Marini",\
                      "shapes": ["triangle", "tetrahedron"],\
                      "orders": [1, 2, 3]},\
                    {"family": "Brezzi-Douglas-Fortin-Marini",\
                      "shapes": ["triangle"],\
                      "orders": [2]},\
                    {"family": "Nedelec 1st kind H(curl)",\
                      "shapes": ["triangle", "tetrahedron"],\
                      "orders": [1, 2, 3]},\
                    {"family": "Nedelec 2nd kind H(curl)",\
                      "shapes": ["triangle", "tetrahedron"],\
                      "orders": [1, 2, 3]},\
                    {"family": "Regge",\
                      "shapes": ["triangle", "tetrahedron"],\
                      "orders": [0, 1, 2, 3]}]

# Create some mixed elements
dg0_tri = FiniteElement("DG", "triangle", 0)
dg1_tri = FiniteElement("DG", "triangle", 1)
cg1_tri = FiniteElement("CG", "triangle", 1)
cr1_tri = FiniteElement("CR", "triangle", 1)
rt1_tri = FiniteElement("RT", "triangle", 1)
drt2_tri = FiniteElement("DRT", "triangle", 2)
bdm1_tri = FiniteElement("BDM", "triangle", 1)
ned1_tri = FiniteElement("N1curl", "triangle", 1)
reg0_tri = FiniteElement("Regge", "triangle", 0)

dg0_tet = FiniteElement("DG", "tetrahedron", 0)
dg1_tet = FiniteElement("DG", "tetrahedron", 1)
cg1_tet = FiniteElement("CG", "tetrahedron", 1)
cr1_tet = FiniteElement("CR", "tetrahedron", 1)
rt1_tet = FiniteElement("RT", "tetrahedron", 1)
drt2_tet = FiniteElement("DRT", "tetrahedron", 2)
bdm1_tet = FiniteElement("BDM", "tetrahedron", 1)
ned1_tet = FiniteElement("N1curl", "tetrahedron", 1)
reg0_tet = FiniteElement("Regge", "tetrahedron", 0)

# Create combinations in pairs.
mix_tri = [MixedElement(e) for e in xcomb([dg0_tri, dg1_tri, cg1_tri, cr1_tri, rt1_tri, drt2_tri, bdm1_tri, ned1_tri, reg0_tri], 2)]
mix_tet = [MixedElement(e) for e in xcomb([dg0_tet, dg1_tet, cg1_tet, cr1_tet, rt1_tet, drt2_tet, bdm1_tet, ned1_tet, reg0_tet], 2)]

mixed_elements = [MixedElement([dg0_tri]*4), MixedElement([cg1_tri]*3), MixedElement([bdm1_tri]*2),\
                  MixedElement([dg1_tri, cg1_tri, cr1_tri, rt1_tri, bdm1_tri, ned1_tri]),\
                  MixedElement([MixedElement([rt1_tri, cr1_tri]), cg1_tri, ned1_tri]),\
                  MixedElement([ned1_tri, dg1_tri, MixedElement([rt1_tri, cr1_tri])]),\
                  MixedElement([drt2_tri, cg1_tri]),\
                  MixedElement([dg0_tet]*4), MixedElement([cg1_tet]*3), MixedElement([bdm1_tet]*2),\
                  MixedElement([dg1_tet, cg1_tet, cr1_tet, rt1_tet, bdm1_tet, ned1_tet]),\
                  MixedElement([MixedElement([rt1_tet, cr1_tet]), cg1_tet, ned1_tet]),\
                  MixedElement([ned1_tet, dg1_tet, MixedElement([rt1_tet, cr1_tet])]),\
                  MixedElement([drt2_tet, cg1_tet]),\
                  MixedElement([cg1_tet, cg1_tet, cg1_tet, reg0_tet])] + mix_tri + mix_tet
