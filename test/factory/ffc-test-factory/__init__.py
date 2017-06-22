# -*- coding: utf-8 -*-

# Copyright (C) 2017 Garth N. Wells
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

import hashlib
import os
import pickle

import ffc_test_factory.factory


def load():

    # Data file
    filename = "element_factory_data.p"
    p = os.path.dirname(os.path.realpath(__file__))
    path = os.path.join(p, filename)

    # Load data file
    elements = pickle.load(open(path, "rb" ))

    # Check hash
    if elements[0] != ffc_test_factory.factory.pickle_hash():
        print("WARNING: Binary code and Python object for factory objects do not match")

    return elements[1]
