"""
This is the compiler, acting as the main interface for compilation
of forms and breaking the compilation into several sequential stages.
The output of each stage is the input of the next stage.

Compiler stage 0: Language, parsing
-----------------------------------

  Input:  Python code or .ufl file
  Output: UFL form

  This stage consists of parsing and expressing a form in the
  UFL form language.

  This stage is completely handled by UFL.

Compiler stage 1: Analysis
--------------------------

  Input:  UFL form
  Output: Preprocessed UFL form and FormData (metadata)

  This stage preprocesses the UFL form and extracts form metadata.
  It may also perform simplifications on the form.

Compiler stage 2: Code representation
-------------------------------------

  Input:  Preprocessed UFL form and FormData (metadata)
  Output: Intermediate Representation (IR)

  This stage examines the input and generates all data needed for code
  generation. This includes generation of finite element basis
  functions, extraction of data for mapping of degrees of freedom and
  possible precomputation of integrals.

  Most of the complexity of compilation is handled in this stage.

  The IR is stored as a dictionary, mapping names of UFC functions to
  data needed for generation of the corresponding code.

Compiler stage 3: Optimization
------------------------------

  Input:  Intermediate Representation (IR)
  Output: Optimized Intermediate Representation (OIR)

  This stage examines the IR and performs optimizations.

  Optimization is currently disabled as a separate stage
  but is implemented as part of the code generation for
  quadrature representation.

Compiler stage 4: Code generation
---------------------------------

  Input:  Optimized Intermediate Representation (OIR)
  Output: C++ code

  This stage examines the OIR and generates the actual C++ code for
  the body of each UFC function.

  The code is stored as a dictionary, mapping names of UFC functions
  to strings containing the C++ code of the body of each function.

Compiler stage 5: Code formatting
---------------------------------

  Input:  C++ code
  Output: C++ code files

  This stage examines the generated C++ code and formats it according
  to the UFC format, generating as output one or more .h/.cpp files
  conforming to the UFC format.

The main interface is defined by the following two functions:

  compile_form
  compile_element

The compiler stages are implemented by the following functions:

  analyze_forms
  or
  analyze_elements  (stage 1)
  compute_ir        (stage 2)
  optimize_ir       (stage 3)
  generate_code     (stage 4)
  format_code       (stage 5)
"""

# Copyright (C) 2007-2013 Anders Logg
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
# Modified by Kristian B. Oelgaard, 2010.
# Modified by Dag Lindbo, 2008.
# Modified by Garth N. Wells, 2009.
# Modified by Martin Alnaes, 2013
#
# First added:  2007-02-05
# Last changed: 2013-01-25

__all__ = ["compile_form", "compile_element"]

# Python modules
from time import time

# FFC modules
from ffc.log import info, info_green, warning
from ffc.parameters import default_parameters

# FFC modules
from ffc.analysis import analyze_forms, analyze_elements
from ffc.representation import compute_ir
from ffc.optimization import optimize_ir
from ffc.codegeneration import generate_code
from ffc.formatting import format_code
from ffc.wrappers import generate_wrapper_code

def compile_form(forms, object_names={}, prefix="Form",\
                 parameters=default_parameters()):
    """This function generates UFC code for a given UFL form or list
    of UFL forms."""

    info("Compiling form %s\n" % prefix)

    # Reset timing
    cpu_time_0 = time()

    # Check input arguments
    forms = _check_forms(forms)
    parameters = _check_parameters(parameters)
    if not forms: return

    # Stage 1: analysis
    cpu_time = time()
    analysis = analyze_forms(forms, object_names, parameters)
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
    code = generate_code(oir, prefix, parameters)
    _print_timing(4, time() - cpu_time)

    # Stage 4.1: generate wrappers
    cpu_time = time()
    wrapper_code = generate_wrapper_code(analysis, prefix, parameters)
    _print_timing(4.1, time() - cpu_time)

    # Stage 5: format code
    cpu_time = time()
    format_code(code, wrapper_code, prefix, parameters)
    _print_timing(5, time() - cpu_time)

    info_green("FFC finished in %g seconds.", time() - cpu_time_0)

def compile_element(elements, prefix="Element", parameters=default_parameters()):
    """This function generates UFC code for a given UFL element or
    list of UFL elements."""

    info("Compiling element %s\n" % prefix)

    # Reset timing
    cpu_time_0 = time()

    # Check input arguments
    elements = _check_elements(elements)
    parameters = _check_parameters(parameters)
    if not elements: return

    # Stage 1: analysis
    cpu_time = time()
    analysis = analyze_elements(elements, parameters)
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
    code = generate_code(oir, prefix, parameters)
    _print_timing(4, time() - cpu_time)

    # Stage 4.1: generate wrappers
    cpu_time = time()
    wrapper_code = generate_wrapper_code(analysis, prefix, parameters)
    _print_timing(4.1, time() - cpu_time)

    # Stage 5: format code
    cpu_time = time()
    format_code(code, wrapper_code, prefix, parameters)
    _print_timing(5, time() - cpu_time)

    info_green("FFC finished in %g seconds.", time() - cpu_time_0)

def _check_forms(forms):
    "Initial check of forms."
    if not isinstance(forms, (list, tuple)):
        forms = (forms,)
    return forms

def _check_elements(elements):
    "Initial check of elements."
    if not isinstance(elements, (list, tuple)):
        elements = (elements,)
    return elements

def _check_parameters(parameters):
    "Initial check of parameters."
    if "blas" in parameters:
        warning("BLAS mode unavailable (will return in a future version).")
    if "quadrature_points" in parameters:
        warning("Option 'quadrature_points' has been replaced by 'quadrature_degree'.")
    return parameters

def _print_timing(stage, timing):
    "Print timing results."
    info("Compiler stage %s finished in %g seconds.\n" % (str(stage), timing))
