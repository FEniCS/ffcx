"This modules provides additional functionality for users of FFC."

# Copyright (C) 2010 Anders Logg
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
# First added:  2010-09-03
# Last changed: 2010-09-06

__all__ = ["compute_tensor_representation"]

# Python modules
from time import time

# FFC modules
from ffc.compiler import _print_timing
from ffc.parameters import default_parameters
from ffc.analysis import analyze_forms
from ffc.representation import compute_ir
from ffc.optimization import optimize_ir
from ffc.codegeneration import generate_code

def compute_tensor_representation(form):
    """Compute tensor representation for given form. This function may
    be useful for those (Hi Matt!) that want to access the FFC tensor
    representation from outside FFC."""

    # Set parameters
    parameters = default_parameters()
    parameters["representation"] = "tensor"
    #parameters["optimize"] = "optimize"

    # The below steps basically duplicate the compiler process but
    # skip the code formatting step. Instead, we extract the relevant
    # portions for tabulate_tensor.

    # Stage 1: analysis
    cpu_time = time()
    analysis = analyze_forms([form], {}, parameters)
    _print_timing(1, time() - cpu_time)

    # Stage 2: intermediate representation
    cpu_time = time()
    ir = compute_ir(analysis, parameters)
    _print_timing(2, time() - cpu_time)

    # Stage 3: optimization
    cpu_time = time()
    oir = optimize_ir(ir, parameters)
    _print_timing(3, time() - cpu_time)

    # Stage 4: code generation
    cpu_time = time()
    code = generate_code(oir, "foo", parameters)
    _print_timing(4, time() - cpu_time)

    # Extract representations
    ir_elements, ir_dofmaps, ir_integrals, ir_forms = ir

    # Extract entries in reference tensor
    reference_tensors = []
    for i in ir_integrals:
        if i["domain_type"] == "cell":
            t = [A0.A0 for (A0, GK, dummy) in i["AK"]]
            if len(t) == 1: t = t[0]
        elif i["domain_type"] == "exterior_facet":
            t = [A0.A0 for j in i["AK"] for (A0, GK, dummy) in j]
            if len(t) == 1: t = t[0]
        elif i["domain_type"] == "interior_facet":
            t = [A0.A0 for j in i["AK"] for k in j for (A0, GK, dummy) in k]
            if len(t) == 1: t = t[0]
        else:
            raise RuntimeError, "Unhandled domain type: %s" % str(i["domain_type"])
        reference_tensors.append(t)

    # Extract code
    code_elements, code_dofmaps, code_integrals, code_forms = code

    # Extract code for computing the geometry tensor
    geometry_tensor_codes = [c["tabulate_tensor"].split("// Compute element tensor")[0] for c in code_integrals]

    # Simplify return values when there is just one term
    if len(reference_tensors) == 1:
        reference_tensors = reference_tensors[0]
    if len(geometry_tensor_codes) == 1:
        geometry_tensor_codes = geometry_tensor_codes[0]

    return reference_tensors, geometry_tensor_codes
