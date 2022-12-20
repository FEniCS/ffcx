# Copyright (C) 2007-2020 Anders Logg and Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Main interface for compilation of forms.

Breaks the compilation into several sequential stages.
The output of each stage is the input of the next stage.

Compiler stages
---------------

0. Language, parsing

   - Input:  Python code or .ufl file
   - Output: UFL form

   This stage consists of parsing and expressing a form in the UFL form
   language. This stage is handled by UFL.

1. Analysis

   - Input:  UFL form
   - Output: Preprocessed UFL form and FormData (metadata)

   This stage preprocesses the UFL form and extracts form metadata. It
   may also perform simplifications on the form.

2. Code representation

   - Input:  Preprocessed UFL form and FormData (metadata)
   - Output: Intermediate Representation (IR)

   This stage examines the input and generates all data needed for code
   generation. This includes generation of finite element basis
   functions, extraction of data for mapping of degrees of freedom and
   possible precomputation of integrals. Most of the complexity of
   compilation is handled in this stage.

   The IR is stored as a dictionary, mapping names of UFC functions to
   data needed for generation of the corresponding code.

3. Code generation

   - Input:  Intermediate Representation (IR)
   - Output: C code

   This stage examines the IR and generates the actual C code for the
   body of each UFC function.

   The code is stored as a dictionary, mapping names of UFC functions to
   strings containing the C code of the body of each function.

4. Code formatting

   - Input:  C code
   - Output: C code files

   This stage examines the generated C++ code and formats it according
   to the UFC format, generating as output one or more .h/.c files
   conforming to the UFC format.

"""

import logging
import typing
from time import time

from ffcx.analysis import analyze_ufl_objects
from ffcx.codegeneration.codegeneration import generate_code
from ffcx.formatting import format_code
from ffcx.ir.representation import compute_ir

logger = logging.getLogger("ffcx")


def _print_timing(stage: int, timing: float):
    logger.info(f"Compiler stage {stage} finished in {timing:.4f} seconds.")


def compile_ufl_objects(ufl_objects: typing.List[typing.Any],
                        object_names: typing.Dict = {},
                        prefix: typing.Optional[str] = None,
                        options: typing.Dict = {},
                        visualise: bool = False):
    """Generate UFC code for a given UFL objects.

    Options
    ----------
    @param ufl_objects:
        Objects to be compiled. Accepts elements, forms, integrals or coordinate mappings.

    """
    # Stage 1: analysis
    cpu_time = time()
    analysis = analyze_ufl_objects(ufl_objects, options)
    _print_timing(1, time() - cpu_time)

    # Stage 2: intermediate representation
    cpu_time = time()
    ir = compute_ir(analysis, object_names, prefix, options, visualise)
    _print_timing(2, time() - cpu_time)

    # Stage 3: code generation
    cpu_time = time()
    code = generate_code(ir, options)
    _print_timing(3, time() - cpu_time)

    # Stage 4: format code
    cpu_time = time()
    code_h, code_c = format_code(code, options)
    _print_timing(4, time() - cpu_time)

    return code_h, code_c
