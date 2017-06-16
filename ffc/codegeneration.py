# -*- coding: utf-8 -*-
"""
Compiler stage 4: Code generation
---------------------------------

This module implements the generation of C++ code for the body of each
UFC function from an (optimized) intermediate representation (OIR).
"""

# Copyright (C) 2009-2017 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, and others
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

from itertools import chain
from ufl import product

# FFC modules
from ffc.log import info, begin, end, debug_code, dstr

# FFC code generation modules
from ffc.representation import pick_representation, ufc_integral_types

import ffc.uflacs.language.cnodes as L
from ffc.uflacs.language.format_lines import format_indented_lines
from ffc.uflacs.backends.ufc.utils import generate_error
from ffc.uflacs.backends.ufc.generators import ufc_integral, ufc_finite_element, ufc_dofmap, ufc_coordinate_mapping, ufc_form

def generate_code(ir, parameters):
    "Generate code from intermediate representation."

    begin("Compiler stage 4: Generating code")

    full_ir = ir

    # FIXME: This has global side effects
    # Set code generation parameters
    # set_float_formatting(parameters["precision"])
    # set_exception_handling(parameters["convert_exceptions_to_warnings"])

    # Extract representations
    ir_finite_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms = ir

    # Generate code for finite_elements
    info("Generating code for %d finite_element(s)" % len(ir_finite_elements))
    code_finite_elements = [_generate_finite_element_code(ir, parameters)
                            for ir in ir_finite_elements]

    # Generate code for dofmaps
    info("Generating code for %d dofmap(s)" % len(ir_dofmaps))
    code_dofmaps = [_generate_dofmap_code(ir, parameters)
                    for ir in ir_dofmaps]

    # Generate code for coordinate_mappings
    info("Generating code for %d coordinate_mapping(s)" % len(ir_coordinate_mappings))
    code_coordinate_mappings = [_generate_coordinate_mapping_code(ir, parameters)
                                for ir in ir_coordinate_mappings]

    # Generate code for integrals
    info("Generating code for integrals")
    code_integrals = [_generate_integral_code(ir, parameters)
                      for ir in ir_integrals]

    # Generate code for forms
    info("Generating code for forms")
    code_forms = [_generate_form_code(ir, parameters)
                  for ir in ir_forms]

    # Extract additional includes
    includes = _extract_includes(full_ir, code_integrals)

    end()

    return (code_finite_elements, code_dofmaps, code_coordinate_mappings,
            code_integrals, code_forms, includes)


def _extract_includes(full_ir, code_integrals):
    ir_finite_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms = full_ir

    # Includes added by representations
    includes = set()
    for code in code_integrals:
        includes.update(code["additional_includes_set"])

    # Includes for dependencies in jit mode
    jit = any(full_ir[i][j]["jit"]
              for i in range(len(full_ir))
              for j in range(len(full_ir[i])))
    if jit:
        dep_includes = set()
        for ir in ir_finite_elements:
            dep_includes.update(_finite_element_jit_includes(ir))
        for ir in ir_dofmaps:
            dep_includes.update(_dofmap_jit_includes(ir))
        for ir in ir_coordinate_mappings:
            dep_includes.update(_coordinate_mapping_jit_includes(ir))
        #for ir in ir_integrals:
        #    dep_includes.update(_integral_jit_includes(ir))
        for ir in ir_forms:
            dep_includes.update(_form_jit_includes(ir))
        includes.update(['#include "%s"' % inc for inc in dep_includes])

    return includes


def _finite_element_jit_includes(ir):
    classnames = ir["create_sub_element"]
    postfix = "_finite_element"
    return [classname.rpartition(postfix)[0] + ".h"
            for classname in classnames]


def _dofmap_jit_includes(ir):
    classnames = ir["create_sub_dofmap"]
    postfix = "_dofmap"
    return [classname.rpartition(postfix)[0] + ".h"
            for classname in classnames]


def _coordinate_mapping_jit_includes(ir):
    classnames = [
        ir["coordinate_finite_element_classname"],
        ir["scalar_coordinate_finite_element_classname"]
        ]
    postfix = "_finite_element"
    return [classname.rpartition(postfix)[0] + ".h"
            for classname in classnames]


def _form_jit_includes(ir):
    # Gather all header names for classes that are separately compiled
    # For finite_element and dofmap the module and header name is the prefix,
    # extracted here with .split, and equal for both classes so we skip dofmap here:
    classnames = list(chain(
        ir["create_finite_element"],
        ir["create_coordinate_finite_element"]
        ))
    postfix = "_finite_element"
    includes = [classname.rpartition(postfix)[0] + ".h"
                for classname in classnames]

    classnames = ir["create_coordinate_mapping"]
    postfix = "_coordinate_mapping"
    includes += [classname.rpartition(postfix)[0] + ".h"
                 for classname in classnames]
    return includes


tt_timing_template = """
    // Initialize timing variables
    static const std::size_t _tperiod = 10000;
    static std::size_t _tcount = 0;
    static auto _tsum = std::chrono::nanoseconds::zero();
    static auto _tavg_best = std::chrono::nanoseconds::max();
    static auto _tmin = std::chrono::nanoseconds::max();
    static auto _tmax = std::chrono::nanoseconds::min();

    // Measure single kernel time
    auto _before = std::chrono::high_resolution_clock::now();
    { // Begin original kernel
%s
    } // End original kernel
    // Measure single kernel time
    auto _after = std::chrono::high_resolution_clock::now();

    // Update time stats
    const std::chrono::seconds _s(1);
    auto _tsingle = _after - _before;
    ++_tcount;
    _tsum += _tsingle;
    _tmin = std::min(_tmin, _tsingle);
    _tmax = std::max(_tmax, _tsingle);

    if (_tcount %% _tperiod == 0 || _tsum > _s)
    {
        // Record best average across batches
        std::chrono::nanoseconds _tavg = _tsum / _tcount;
        if (_tavg_best > _tavg)
            _tavg_best = _tavg;

        // Convert to ns
        auto _tot_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(_tsum).count();
        auto _avg_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(_tavg).count();
        auto _min_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(_tmin).count();
        auto _max_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(_tmax).count();
        auto _avg_best_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(_tavg_best).count();

        // Print report
        std::cout << "FFC tt time:"
                  << "  avg_best = " << _avg_best_ns << " ns,"
                  << "  avg = " << _avg_ns << " ns,"
                  << "  min = " << _min_ns << " ns,"
                  << "  max = " << _max_ns << " ns,"
                  << "  tot = " << _tot_ns << " ns,"
                  << "  n = " << _tcount
                  << std::endl;

        // Reset statistics for next batch
        _tcount = 0;
        _tsum = std::chrono::nanoseconds(0);
        _tmin = std::chrono::nanoseconds::max();
        _tmax = std::chrono::nanoseconds::min();
    }
"""


def _generate_tabulate_tensor_comment(ir, parameters):
    "Generate comment for tabulate_tensor."

    r = ir["representation"]
    integrals_metadata = ir["integrals_metadata"]
    integral_metadata = ir["integral_metadata"]

    comment  = [ L.Comment("This function was generated using '%s' representation" % r),
                 L.Comment("with the following integrals metadata:"),
                 L.Comment(""),
                 L.Comment("\n".join(dstr(integrals_metadata).split("\n")[:-1]))
                 ]
    for i, metadata in enumerate(integral_metadata):
        comment += [ L.Comment(""),
                     L.Comment("and the following integral %d metadata:" % i),
                     L.Comment(""),
                     L.Comment("\n".join(dstr(metadata).split("\n")[:-1]))
                   ]

    return format_indented_lines(L.StatementList(comment).cs_format(0), 1)

def _generate_integral_code(ir, parameters):
    "Generate code for integrals from intermediate representation."
    # Select representation
    r = pick_representation(ir["representation"])

    # Generate code
    # TODO: Drop prefix argument and get from ir:
    code = r.generate_integral_code(ir, ir["prefix"], parameters)

    # Hack for benchmarking overhead in assembler with empty
    # tabulate_tensor
    if parameters["generate_dummy_tabulate_tensor"]:
        code["tabulate_tensor"] = ""

    # Wrapping tabulate_tensor in a timing snippet for benchmarking
    if parameters["add_tabulate_tensor_timing"]:
        code["tabulate_tensor"] = tt_timing_template % code["tabulate_tensor"]
        code["additional_includes_set"] = code.get("additional_includes_set", set())
        code["additional_includes_set"].add("#include <chrono>")
        code["additional_includes_set"].add("#include <iostream>")

    # Generate comment
    code["tabulate_tensor_comment"] = _generate_tabulate_tensor_comment(ir, parameters)

    return code

# TODO: Replace the above with this, currently something is not working
def _new_generate_integral_code(ir, parameters):
    "Generate code for integrals from intermediate representation."
    return ufc_integral(ir["integral_type"]).generate_snippets(L, ir, parameters)


def _generate_finite_element_code(ir, parameters):
    "Generate code for finite_element from intermediate representation."
    return ufc_finite_element().generate_snippets(L, ir, parameters)


def _generate_dofmap_code(ir, parameters):
    "Generate code for dofmap from intermediate representation."
    return ufc_dofmap().generate_snippets(L, ir, parameters)


def _generate_coordinate_mapping_code(ir, parameters):
    "Generate code for coordinate_mapping from intermediate representation."
    return ufc_coordinate_mapping().generate_snippets(L, ir, parameters)


def _generate_form_code(ir, parameters):
    "Generate code for coordinate_mapping from intermediate representation."
    return ufc_form().generate_snippets(L, ir, parameters)
