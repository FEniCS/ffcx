# -*- coding: utf-8 -*-
# Copyright (C) 2007-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Main interface for compilation
of forms and breaking the compilation into several sequential stages.
The output of each stage is the input of the next stage.

Compiler stages:

#. Language, parsing

   - Input:  Python code or .ufl file
   - Output: UFL form

   This stage consists of parsing and expressing a form in the UFL form
   language. This stage is handled by UFL.

#. Analysis

   - Input:  UFL form
   - Output: Preprocessed UFL form and FormData (metadata)

   This stage preprocesses the UFL form and extracts form metadata. It
   may also perform simplifications on the form.

#. Code representation

   - Input:  Preprocessed UFL form and FormData (metadata)
   - Output: Intermediate Representation (IR)

   This stage examines the input and generates all data needed for code
   generation. This includes generation of finite element basis
   functions, extraction of data for mapping of degrees of freedom and
   possible precomputation of integrals. Most of the complexity of
   compilation is handled in this stage.

   The IR is stored as a dictionary, mapping names of UFC functions to
   data needed for generation of the corresponding code.

#. Code generation

   - Input:  Intermediate Representation (IR)
   - Output: C code

   This stage examines the IR and generates the actual C code for the
   body of each UFC function.

   The code is stored as a dictionary, mapping names of UFC functions to
   strings containing the C code of the body of each function.

#. Code formatting

   - Input:  C code
   - Output: C code files

   This stage examines the generated C++ code and formats it according
   to the UFC format, generating as output one or more .h/.c files
   conforming to the UFC format.

"""

import logging
import os
from collections import defaultdict
from time import time
from typing import Dict, List, Tuple, Union

from ffc.analysis import analyze_ufl_objects
from ffc.codegeneration.codegeneration import generate_code
from ffc.formatting import format_code
from ffc.parameters import validate_parameters
from ffc.ir.representation import compute_ir
from ffc.wrappers import generate_wrapper_code

logger = logging.getLogger(__name__)


def _print_timing(stage, timing):
    logger.info("Compiler stage {stage} finished in {time} seconds.".format(
        stage=stage, time=timing))


def compile_ufl_objects(ufl_objects: Union[List, Tuple],
                        object_names: Dict = {},
                        prefix: str = None,
                        parameters: Dict = None,
                        jit: bool = False):
    """Generate UFC code for a given UFL objects.

    Parameters
    ----------
    ufl_objects
        Objects to be compiled. Accepts elements, forms, integrals or coordinate mappings.

    """
    logger.info("Compiling {}\n".format(prefix))

    # Reset timing
    cpu_time_0 = time()

    # Note that jit will always pass validated parameters so
    # this is only for commandline and direct call from python
    if not jit:
        parameters = validate_parameters(parameters)

    # Check input arguments
    if not isinstance(ufl_objects, (list, tuple)):
        ufl_objects = (ufl_objects, )
    if not ufl_objects:
        return "", ""

    if prefix != os.path.basename(prefix):
        raise RuntimeError("Invalid prefix, looks like a full path? prefix='{}'.".format(prefix))

    # Check that all UFL objects passed here are of the same class/type
    obj_type = type(ufl_objects[0])
    assert (obj_type == type(x) for x in ufl_objects)

    # Stage 1: analysis
    cpu_time = time()
    analysis = analyze_ufl_objects(ufl_objects, parameters)
    _print_timing(1, time() - cpu_time)

    # Stage 2: intermediate representation
    cpu_time = time()
    ir = compute_ir(analysis, prefix, parameters, jit)
    _print_timing(2, time() - cpu_time)

    # Stage 4: code generation
    cpu_time = time()
    code = generate_code(ir, parameters, jit)
    _print_timing(4, time() - cpu_time)

    # Stage 4.1: generate convenience wrappers, e.g. for DOLFIN
    cpu_time = time()

    # Extract class names from the IR and add to a dict
    # ir_finite_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms = ir
    if len(object_names) > 0:
        classnames = defaultdict(list)
        comp = ["elements", "dofmaps", "coordinate_maps", "integrals", "forms"]
        for ir_comp, e_name in zip(ir, comp):
            for e in ir_comp:
                classnames[e_name].append(e["classname"])
        wrapper_code = generate_wrapper_code(analysis, prefix, object_names, classnames, parameters)
    else:
        wrapper_code = None

    _print_timing(4.1, time() - cpu_time)

    # Stage 5: format code
    cpu_time = time()
    code_h, code_c = format_code(code, wrapper_code, prefix, parameters)
    _print_timing(5, time() - cpu_time)

    logger.info("FFC finished in {} seconds.".format(time() - cpu_time_0))

    return code_h, code_c
