# Copyright (C) 2013 Martin S. Alnaes
"""Cell geometry demo.

A functional M involving a bunch of cell geometry quantities.
"""

import basix.ufl
from ufl import (
    CellVolume,
    Circumradius,
    Coefficient,
    FacetArea,
    FacetNormal,
    FunctionSpace,
    Mesh,
    SpatialCoordinate,
    TrialFunction,
    ds,
    dx,
)
from ufl.geometry import FacetEdgeVectors

V = basix.ufl.element("P", "tetrahedron", 1)
domain = Mesh(basix.ufl.element("P", "tetrahedron", 1, shape=(3,)))
space = FunctionSpace(domain, V)
u = Coefficient(space)

# TODO: Add all geometry for all cell types to this and other demo
# files, need for regression test.
x = SpatialCoordinate(domain)
n = FacetNormal(domain)
vol = CellVolume(domain)
rad = Circumradius(domain)
area = FacetArea(domain)

M = u * (x[0] * vol * rad) * dx + u * (x[0] * vol * rad * area) * ds

# Test some obscure functionality
fev = FacetEdgeVectors(domain)
v = TrialFunction(space)
L = fev[0, 0] * v * ds
