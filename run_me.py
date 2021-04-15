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

import basix
basix.create_element("Nedelec 1st kind H(curl)", "tetrahedron", 1)


print("A")
#element = create_basix_element(FiniteElement("N1curl", "tetrahedron", 1))
print("B")
