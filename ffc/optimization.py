# -*- coding: utf-8 -*-

# Copyright (C) 2009-2017 Anders Logg
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
# Modified by Martin Sandve Alnæs, 2013

"""
Compiler stage 5: optimization
------------------------------

This module implements the optimization of an intermediate code
representation.
"""

# FFC modules
from ffc.log import info, begin, end
from ffc.representation import pick_representation


def optimize_ir(ir, parameters):
    "Optimize intermediate form representation."

    begin("Compiler stage 3: Optimizing intermediate representation")

    # Extract representations
    ir_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms = ir

    # Check if optimization is requested
    if not any(ir["integrals_metadata"]["optimize"] for ir in ir_integrals):
        info(r"Skipping optimizations, add -O or attach {'optimize': True} "
             "metadata to integrals")

    # Call on every bunch of integrals wich are compiled together
    oir_integrals = [_optimize_integral_ir(ir, parameters)
                     if ir["integrals_metadata"]["optimize"] else ir
                     for ir in ir_integrals]

    end()

    return ir_elements, ir_dofmaps, ir_coordinate_mappings, oir_integrals, ir_forms


def _optimize_integral_ir(ir, parameters):
    "Compute optimized intermediate represention of integral."
    r = pick_representation(ir["representation"])
    return r.optimize_integral_ir(ir, parameters)
