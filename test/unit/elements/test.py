"""Unit tests for FFC finite elements"""

# Copyright (C) 2013 Marie E. Rognes
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

import unittest

from ufl import interval
from ufl import FiniteElement

from ffc import compile_element

class TestCompileElements(unittest.TestCase):

    def testRadau(self):
        "Test that Radau elements compile."
        for degree in range(3):
            element = FiniteElement("Radau", interval, degree)
            compile_element(element)

    def testLobatto(self):
        "Test that Lobatto elements compile."
        for degree in range(1, 4):
            element = FiniteElement("Lobatto", interval, degree)
            compile_element(element)


if __name__ == "__main__":
    unittest.main()
