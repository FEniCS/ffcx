# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# import tempfile
# from pathlib import Path
# import importlib
import logging
import os
import subprocess
# import sys
# import time
import re

# import cffi
import llvmlite.binding as llvm

import ffc
import ffc.config

logger = logging.getLogger(__name__)

UFC_HEADER_DECL = """
typedef {} ufc_scalar_t;  /* Hack to deal with scalar type */

typedef struct ufc_coordinate_mapping ufc_coordinate_mapping;
typedef struct ufc_finite_element ufc_finite_element;
typedef struct ufc_dofmap ufc_dofmap;

typedef enum
{{
interval = 10,
triangle = 20,
quadrilateral = 30,
tetrahedron = 40,
hexahedron = 50,
vertex = 60,
}} ufc_shape;
"""

# Get declarations directly from ufc.h
file_dir = os.path.dirname(os.path.abspath(__file__))
with open(file_dir + "/ufc.h", "r") as f:
    ufc_h = ''.join(f.readlines())

UFC_ELEMENT_DECL = '\n'.join(re.findall('typedef struct ufc_finite_element.*?ufc_finite_element;', ufc_h, re.DOTALL))
UFC_DOFMAP_DECL = '\n'.join(re.findall('typedef struct ufc_dofmap.*?ufc_dofmap;', ufc_h, re.DOTALL))
UFC_COORDINATEMAPPING_DECL = '\n'.join(re.findall('typedef struct ufc_coordinate_mapping.*?ufc_coordinate_mapping;',
                                                  ufc_h, re.DOTALL))
UFC_FORM_DECL = '\n'.join(re.findall('typedef struct ufc_form.*?ufc_form;', ufc_h, re.DOTALL))
UFC_INTEGRAL_DECL = '\n'.join(re.findall('typedef struct ufc_integral.*?ufc_integral;', ufc_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall('typedef struct ufc_custom_integral.*?ufc_custom_integral;',
                                          ufc_h, re.DOTALL))


def compile_elements(elements, parameters=None):
    """Compile a list of UFL elements and dofmaps into UFC Python objects"""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling elements: ' + str(elements))

    # Get a signature for these elements
    module_name = 'libffc_elements_' + ffc.classname.compute_signature(elements, '', p)

    names = []
    for e in elements:
        name = ffc.ir.representation.make_finite_element_jit_classname(e, "JIT", p)
        names.append(name)
        name = ffc.ir.representation.make_dofmap_jit_classname(e, "JIT", p)
        names.append(name)

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL
    element_template = "ufc_finite_element * create_{name}(void);\n"
    dofmap_template = "ufc_dofmap * create_{name}(void);\n"

    for i in range(len(elements)):
        decl += element_template.format(name=names[i * 2])
        decl += dofmap_template.format(name=names[i * 2 + 1])

    objects, module = _compile_objects(decl, elements, names, module_name, p)
    # Pair up elements with dofmaps
    objects = list(zip(objects[::2], objects[1::2]))
    return objects, module


def compile_forms(forms, parameters=None):
    """Compile a list of UFL forms into UFC Python objects"""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling forms: ' + str(forms))

    # Get a signature for these forms
    module_name = 'libffc_forms_' + ffc.classname.compute_signature(forms, '', p)

    form_names = [ffc.classname.make_name("JIT", "form", i)
                  for i in range(len(forms))]

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL + \
        UFC_COORDINATEMAPPING_DECL + UFC_INTEGRAL_DECL + UFC_FORM_DECL

    form_template = "ufc_form * create_{name}(void);\n"
    for name in form_names:
        decl += form_template.format(name=name)

    return _compile_objects(decl, forms, form_names, module_name, p)


def compile_coordinate_maps(meshes, parameters=None):
    """Compile a list of UFL coordinate mappings into UFC Python objects"""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling cmaps: ' + str(meshes))

    # Get a signature for these cmaps
    module_name = 'libffc_cmaps_' + ffc.classname.compute_signature(meshes, '', p, True)

    cmap_names = [ffc.ir.representation.make_coordinate_mapping_jit_classname(
        mesh.ufl_coordinate_element(), "JIT", p) for mesh in meshes]

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_COORDINATEMAPPING_DECL
    cmap_template = "ufc_coordinate_mapping * create_{name}(void);\n"

    for name in cmap_names:
        decl += cmap_template.format(name=name)

    return _compile_objects(decl, meshes, cmap_names, module_name, p)


def _compile_objects(decl, ufl_objects, object_names, module_name, parameters):

    _, code_body = ffc.compiler.compile_ufl_objects(ufl_objects, prefix="JIT", parameters=parameters)

    command = "clang -x c - {includes} -c -S -emit-llvm -o -".format(includes="-I"
                                                                     + ffc.codegeneration.get_include_path())

    ps = subprocess.Popen(command.split(" "), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    code_body = code_body.encode('utf-8')
    llvm_ir = ps.communicate(input=code_body)
    llvm_ir = llvm_ir[0].decode('utf-8')

    #    print(llvm_ir)

    # All these initializations are required for code generation!
    llvm.initialize()
    llvm.initialize_native_target()
    llvm.initialize_native_asmprinter()

    # Create a target machine representing the host
    target = llvm.Target.from_default_triple()
    target_machine = target.create_target_machine()
    # And an execution engine with an empty backing module
    backing_mod = llvm.parse_assembly("")
    engine = llvm.create_mcjit_compiler(backing_mod, target_machine)

    mod = llvm.parse_assembly(llvm_ir)
    mod.verify()
    # Now add the module and make sure it is ready for execution
    engine.add_module(mod)
    engine.finalize_object()
    engine.run_static_constructors()

    fnames = ['create_' + name for name in object_names]

    for f in fnames:
        q = engine.get_function_address(f)
        print(f, q)

    quit()

    # Build list of compiled objects
    #    compiled_objects = [getattr(compiled_module.lib, "create_" + name)() for name in object_names]

    compiled_objects = None
    compiled_module = None

    return compiled_objects, compiled_module
