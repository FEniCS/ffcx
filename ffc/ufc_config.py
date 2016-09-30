# Copyright (C) 2016 Jan Blechta, Johannes Ring
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


# TODO: Can we safely remove this?
def get_ufc_cxx_flags():
    """Return C++ flags for compiling UFC C++11 code. Return type
    is a list of strings.
    """
    return ["-std=c++11"]


# TODO: Can we safely remove this name?
from ffc.backends.ufc import get_ufc_signature


# ufc_signature() already introduced to FFC standard in 1.7.0dev,
# called by the dolfin cmake build system to compare against
# future imported ffc versions for compatibility.
from ffc.backends.ufc import get_ufc_signature as ufc_signature
