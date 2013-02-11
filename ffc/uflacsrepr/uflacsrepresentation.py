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
# Last changed: 2013-02-11

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
    ffc_data = None

    # Call upon flacs to build graphs and ssa or other intermediate representation and add to ir
    from uflacs.backends.ffc import compute_uflacs_integral_ir
    uir = compute_uflacs_integral_ir(ffc_data, itg_data, form_data, parameters)
    ir.update(uir)

    return ir

def optimize_integral_ir(ir, parameters):
    "Compute optimized intermediate representation of integral."

    info("Optimizing uflacs representation")

    # Call upon uflacs to optimize ssa representation prior to code generation. Should be possible to skip this step.
    from uflacs.backends.ffc import optimize_integral_ir
    oir = optimize_integral_ir(ir, parameters)

    return oir

def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from uflacs representation")

    # Generate code
    code = initialize_integral_code(ir, prefix, parameters)

    # Call upon uflacs to render graphs into code
    from uflacs.backends.ffc import generate_tabulate_tensor_code
    code["tabulate_tensor"] = generate_tabulate_tensor_code(ir, parameters)

    code["tabulate_tensor_quadrature"] = format["do nothing"] # TODO: Remove
    return code
