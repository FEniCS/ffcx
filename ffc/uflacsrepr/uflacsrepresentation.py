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
# Last changed: 2013-02-12

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

    # TODO: Call upon ffc to build ir for element tables etc
    ffc_data = None

    # Delegate to flacs to build its intermediate representation and add to ir
    import uflacs.backends.ffc
    uir = uflacs.backends.ffc.compute_tabulate_tensor_ir(ffc_data, itg_data, form_data, parameters)
    ir.update(uir)

    return ir

def optimize_integral_ir(ir, parameters):
    "Compute optimized intermediate representation of integral."

    info("Optimizing uflacs representation")

    # Call upon uflacs to optimize ssa representation prior to code generation. Should be possible to skip this step.
    import uflacs.backends.ffc
    oir = uflacs.backends.ffc.optimize_tabulate_tensor_ir(ir, parameters)

    return oir

def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from uflacs representation")

    # Generate generic ffc code snippets
    code = initialize_integral_code(ir, prefix, parameters)

    # TODO: Call upon ffc to generate code for element tables etc
    #code["..."] = ...

    # Delegate to uflacs to generate tabulate_tensor body
    import uflacs.backends.ffc
    code["tabulate_tensor"] = uflacs.backends.ffc.generate_tabulate_tensor_code(ir, parameters)

    code["tabulate_tensor_quadrature"] = format["do nothing"] # TODO: Remove
    return code
