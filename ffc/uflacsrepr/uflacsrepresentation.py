# Copyright (C) 2013 Martin Alnaes
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
# First added:  2009-12-16
# Last changed: 2013-02-10

from ffc.representationutils import initialize_integral_ir, initialize_integral_code
from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning
from ffc.cpp import format

def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing uflacs representation")

    # Initialise representation
    ir = initialize_integral_ir("uflacs", itg_data, form_data, form_id)

    # TODO: Call upon ffc to build element tabules or other intermediate representation and add to ir
    #ffc_data = ...

    # TODO: Call upon flacs to build graphs and ssa or other intermediate representation and add to ir
    #import uflacs
    uir = {} #uflacs.backends.ffc.compute_integral_ir(itg_data, form_data, ffc_data, parameters)
    ir.update(uir)

    # FIXME: Return something minimal that will pass through ffc
    return ir

def optimize_integral_ir(ir, parameters):
    "Compute optimized intermediate representation of integral."

    info("Optimizing uflacs representation")

    # TODO: Call upon uflacs to optimize ssa representation prior to code generation. Should be possible to skip this step.

    return ir

def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from uflacs representation")

    # Generate code
    code = initialize_integral_code(ir, prefix, parameters)

    code["tabulate_tensor"] = format["do nothing"]
    code["tabulate_tensor_quadrature"] = format["do nothing"]

    # TODO: Call upon uflacs to render graphs into code
    #import uflacs
    #code["tabulate_tensor"] = uflacs.backends.ffc.generate_integral_code(ir, parameters)

    # FIXME: Return something minimal that will pass through ffc
    return code
