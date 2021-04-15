# Copyright (C) 2007-2017 Anders Logg and Garth N. Wells
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
# Modified by Marie E. Rognes, 2010
# Modified by Lizao Li, 2016


import numpy

from ffcx.basix_interface import create_basix_element
from ufl import FiniteElement


def element_coords(cell):
    if cell == "interval":
        return [(0,), (1,)]
    elif cell == "triangle":
        return [(0, 0), (1, 0), (0, 1)]
    elif cell == "tetrahedron":
        return [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    elif cell == "quadrilateral":
        return [(0, 0), (1, 0), (0, 1), (1, 1)]
    elif cell == "hexahedron":
        return [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
    else:
        RuntimeError("Unknown cell type")


def random_point(shape):
    w = numpy.random.random(len(shape))
    return sum([numpy.array(shape[i]) * w[i] for i in range(len(shape))]) / sum(w)


print("A")
element = create_basix_element(FiniteElement("N1curl", "tetrahedron", 1))
print("B")
points = [random_point(element_coords(cell)) for i in range(5)]
print("C")
for i, x in enumerate(points):
    print(f"D{i}")
    table = element.tabulate(0, (x,))
    print(f"E{i}")
    basis = table[0]
    print(f"F{i}")
    for j in range(3):
        print(f"G{i},{j}")
        print(basis[0][j::3])
        print(f"H{i},{j}")
