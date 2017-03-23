# -*- coding: utf-8 -*-
"""
Compiler stage 4: Code generation
---------------------------------

This module implements the generation of C++ code for the body of each
UFC function from an (optimized) intermediate representation (OIR).
"""

# Copyright (C) 2009-2016 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
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
from ffc.cpp import format, indent
from ffc.cpp import set_exception_handling, set_float_formatting

# FFC code generation modules
from ffc.evaluatebasis import _evaluate_basis, _evaluate_basis_all
from ffc.evaluatebasisderivatives import _evaluate_basis_derivatives
from ffc.evaluatebasisderivatives import _evaluate_basis_derivatives_all
from ffc.evaluatedof import evaluate_dof_and_dofs, affine_weights
from ffc.interpolatevertexvalues import interpolate_vertex_values
from ffc.representation import pick_representation, ufc_integral_types


# Errors issued for non-implemented functions
def _not_implemented(function_name, return_null=False):
    body = format["exception"]("%s not yet implemented." % function_name)
    if return_null:
        body += "\n" + format["return"](0)
    return body


def generate_code(ir, parameters):
    "Generate code from intermediate representation."

    begin("Compiler stage 4: Generating code")

    # FIXME: This has global side effects
    # Set code generation parameters
    set_float_formatting(parameters["precision"])
    set_exception_handling(parameters["convert_exceptions_to_warnings"])

    # Extract representations
    ir_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms = ir

    # Generate code for elements
    info("Generating code for %d element(s)" % len(ir_elements))
    code_elements = [_generate_element_code(ir, parameters)
                     for ir in ir_elements]

    # Generate code for dofmaps
    info("Generating code for %d dofmap(s)" % len(ir_dofmaps))
    code_dofmaps = [_generate_dofmap_code(ir, parameters)
                    for ir in ir_dofmaps]

    # Generate code for coordinate_mappings
    info("Generating code for %d coordinate_mapping(s)" % len(ir_coordinate_mappings))
    code_coordinate_mappings = [_generate_coordinate_mapping_code(ir, parameters)
                                for ir in ir_coordinate_mappings]
    # FIXME: This disables output of generated coordinate_mapping class, until implemented properly
    code_coordinate_mappings = []

    # Generate code for integrals
    info("Generating code for integrals")
    code_integrals = [_generate_integral_code(ir, parameters)
                      for ir in ir_integrals]

    # Generate code for forms
    info("Generating code for forms")
    code_forms = [_generate_form_code(ir, parameters)
                  for ir in ir_forms]

    end()

    return code_elements, code_dofmaps, code_coordinate_mappings, code_integrals, code_forms


def _generate_element_code(ir, parameters):
    "Generate code for finite element from intermediate representation."

    # Skip code generation if ir is None
    if ir is None:
        return None

    # Prefetch formatting to speedup code generation
    ret = format["return"]
    do_nothing = format["do nothing"]
    create = format["create foo"]

    # Codes generated together
    (evaluate_dof_code, evaluate_dofs_code) \
        = evaluate_dof_and_dofs(ir["evaluate_dof"])

    #element_number = ir["id"]

    # Generate code
    code = {}
    code["classname"] = ir["classname"]
    code["members"] = ""
    code["constructor"] = do_nothing
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = do_nothing

    code["signature"] = ret('"%s"' % ir["signature"])
    code["cell_shape"] = ret(format["cell"](ir["cell_shape"]))
    code["topological_dimension"] = ret(ir["topological_dimension"])
    code["geometric_dimension"] = ret(ir["geometric_dimension"])
    code["space_dimension"] = ret(ir["space_dimension"])

    code["value_rank"] = ret(len(ir["value_shape"]))
    code["value_dimension"] = _value_dimension(ir["value_shape"])
    code["value_size"] = ret(product(ir["value_shape"]))
    code["reference_value_rank"] = ret(len(ir["reference_value_shape"]))
    code["reference_value_dimension"] = _value_dimension(ir["reference_value_shape"])
    code["reference_value_size"] = ret(product(ir["reference_value_shape"]))

    code["degree"] = ret(ir["degree"])
    code["family"] = ret('"%s"' % (ir["family"],))

    code["evaluate_basis"] = _evaluate_basis(ir["evaluate_basis"])
    code["evaluate_basis_all"] = _evaluate_basis_all(ir["evaluate_basis"])
    code["evaluate_basis_derivatives"] \
        = _evaluate_basis_derivatives(ir["evaluate_basis"])
    code["evaluate_basis_derivatives_all"] \
        = _evaluate_basis_derivatives_all(ir["evaluate_basis"])

    code["evaluate_dof"] = evaluate_dof_code
    code["evaluate_dofs"] = evaluate_dofs_code
    code["interpolate_vertex_values"] \
        = interpolate_vertex_values(ir["interpolate_vertex_values"])
    code["tabulate_dof_coordinates"] \
        = _tabulate_dof_coordinates(ir["tabulate_dof_coordinates"])
    code["num_sub_elements"] = ret(ir["num_sub_elements"])
    code["create_sub_element"] = _create_sub_element(ir)
    code["create"] = ret(create(code["classname"]))
    code["additional_includes_set"] = _additional_includes_finite_element(ir)

    # Postprocess code
    _postprocess_code(code, parameters)

    return code


def _generate_dofmap_code(ir, parameters):
    "Generate code for dofmap from intermediate representation."

    # Skip code generation if ir is None
    if ir is None:
        return None

    # Prefetch formatting to speedup code generation
    ret = format["return"]
    declare = format["declaration"]
    assign = format["assign"]
    do_nothing = format["do nothing"]
    switch = format["switch"]
    f_int = format["int"]
    f_d = format["argument dimension"]
    create = format["create foo"]

    #element_number = ir["id"]

    # Generate code
    code = {}
    code["classname"] = ir["classname"]
    code["members"] = ""
    code["constructor"] = do_nothing
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = do_nothing
    code["signature"] = ret('"%s"' % ir["signature"])
    code["needs_mesh_entities"] \
        = _needs_mesh_entities(ir["needs_mesh_entities"])
    code["topological_dimension"] = ret(ir["topological_dimension"])
    code["global_dimension"] = _global_dimension(ir["global_dimension"])
    code["num_element_dofs"] = ret(ir["num_element_dofs"])
    code["num_facet_dofs"] = ret(ir["num_facet_dofs"])
    code["num_entity_dofs"] \
        = switch(f_d, [ret(num) for num in ir["num_entity_dofs"]],
                 ret(f_int(0)))
    code["num_entity_closure_dofs"] \
        = switch(f_d, [ret(num) for num in ir["num_entity_closure_dofs"]],
                 ret(f_int(0)))
    code["tabulate_dofs"] = _tabulate_dofs(ir["tabulate_dofs"])
    code["tabulate_facet_dofs"] \
        = _tabulate_facet_dofs(ir["tabulate_facet_dofs"])
    code["tabulate_entity_dofs"] \
        = _tabulate_entity_dofs(ir["tabulate_entity_dofs"])
    code["tabulate_entity_closure_dofs"] \
        = _tabulate_entity_closure_dofs(ir["tabulate_entity_closure_dofs"])
    code["num_sub_dofmaps"] = ret(ir["num_sub_dofmaps"])
    code["create_sub_dofmap"] = _create_sub_dofmap(ir)
    code["create"] = ret(create(code["classname"]))
    code["additional_includes_set"] = _additional_includes_dofmap(ir)

    # Postprocess code
    _postprocess_code(code, parameters)

    return code


def _additional_includes_dofmap(ir):
    if not ir["jit"]:
        return set()
    dofmap_classnames = ir["create_sub_dofmap"]
    jit_includes = [classname.split("_dofmap")[0] + ".h"
                    for classname in dofmap_classnames]
    return set("#include <%s>" % inc for inc in jit_includes)


def _additional_includes_finite_element(ir):
    if not ir["jit"]:
        return set()
    finite_element_classnames = ir["create_sub_element"]
    jit_includes = [classname.split("_finite_element")[0] + ".h"
                    for classname in finite_element_classnames]
    return set("#include <%s>" % inc for inc in jit_includes)


def _additional_includes_coordinate_mapping(ir):
    if not ir["jit"]:
        return set()
    finite_element_classnames = [
        ir["coordinate_finite_element_classname"],
        ir["scalar_coordinate_finite_element_classname"]
        ]
    jit_includes = [classname.split("_finite_element")[0] + ".h"
                    for classname in finite_element_classnames]
    return set("#include <%s>" % inc for inc in jit_includes)


def _additional_includes_form(ir):
    if not ir["jit"]:
        return set()
    # Gather all header names for classes that are separately compiled
    # For finite_element and dofmap the module and header name is the prefix,
    # extracted here with .split, and equal for both classes so we skip dofmap here:
    finite_element_classnames = list(chain(
        ir["create_finite_element"],
        ir["create_coordinate_finite_element"]
        ))
    jit_includes = set(classname.split("_finite_element")[0] + ".h"
                       for classname in finite_element_classnames)
    # FIXME: Enable when coordinate_mapping is fully generated:
    #jit_includes.update(classname + ".h" for classname in ir["create_coordinate_mapping"])
    return set("#include <%s>" % inc for inc in jit_includes)


def _generate_coordinate_mapping_code(ir, parameters):
    "Generate code for coordinate_mapping from intermediate representation."

    # Skip code generation if ir is None
    if ir is None:
        return None

    coordinate_mapping_number = ir["id"]

    # FIXME: Get code dict from current work in uflacs

    # Generate code
    code = {}
    code["classname"] = ir["classname"]

    code["members"] = ""
    code["constructor"] = ""
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = ""
    code["create"] = "return nullptr;"

    code["signature"] = 'return "fixme";'
    code["cell_shape"] = "return ufc::shape::triangle;"
    code["topological_dimension"] = "return 2;"
    code["geometric_dimension"] = "return 3;"

    code["create_coordinate_finite_element"] = "return nullptr;"
    code["create_coordinate_dofmap"] = "return nullptr;"

    code["compute_physical_coordinates"] = ""
    code["compute_reference_coordinates"] = ""
    code["compute_jacobians"] = ""
    code["compute_jacobian_determinants"] = ""
    code["compute_jacobian_inverses"] = ""
    code["compute_geometry"] = ""

    code["additional_includes_set"] = _additional_includes_coordinate_mapping(ir)

    return code


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


def _generate_integral_code(ir, parameters):
    "Generate code for integrals from intermediate representation."

    # Skip code generation if ir is None
    if ir is None:
        return None

    # Select representation
    r = pick_representation(ir["representation"])

    # Generate code
    # TODO: Drop prefix argument and get from ir:
    code = r.generate_integral_code(ir, ir["prefix"], parameters)

    # Hack for benchmarking overhead in assembler with empty tabulate_tensor
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

    # Indent code (unused variables should already be removed)
    # FIXME: Remove this quick hack
    if ir["representation"] != "uflacs":
        _indent_code(code)
    else:
        code["tabulate_tensor_comment"] = indent(code["tabulate_tensor_comment"], 4)

    return code


def _generate_tabulate_tensor_comment(ir, parameters):
    "Generate comment for tabulate_tensor."

    r = ir["representation"]
    integrals_metadata = ir["integrals_metadata"]
    integral_metadata = ir["integral_metadata"]

    comment  = format["comment"]("This function was generated using '%s' representation" % r) + "\n"
    comment += format["comment"]("with the following integrals metadata:") + "\n"
    comment += format["comment"]("") + "\n"
    comment += "\n".join([format["comment"]("  " + l) for l in dstr(integrals_metadata).split("\n")][:-1])
    comment += "\n"
    for i, metadata in enumerate(integral_metadata):
        comment += format["comment"]("") + "\n"
        comment += format["comment"]("and the following integral %d metadata:" % i) + "\n"
        comment += format["comment"]("") + "\n"
        comment += "\n".join([format["comment"]("  " + l) for l in dstr(metadata).split("\n")][:-1])
        comment += "\n"

    return comment


def _generate_original_coefficient_position(original_coefficient_positions):
    # TODO: I don't know how to implement this using the format dict,
    # this will do for now:
    initializer_list = ", ".join(str(i)
                                 for i in original_coefficient_positions)
    code = '\n'.join([
        "static const std::vector<std::size_t> position({%s});"
        % initializer_list, "return position[i];",
    ])
    return code


def _generate_form_code(ir, parameters):
    "Generate code for form from intermediate representation."

    # Skip code generation if ir is None
    if ir is None:
        return None

    # Prefetch formatting to speedup code generation
    ret = format["return"]
    do_nothing = format["do nothing"]

    # Generate code
    code = {}
    code["classname"] = ir["classname"]
    code["members"] = ""

    code["constructor"] = do_nothing
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = do_nothing

    code["signature"] = ret('"%s"' % ir["signature"])
    code["original_coefficient_position"] = \
        _generate_original_coefficient_position(ir["original_coefficient_position"])
    code["rank"] = ret(ir["rank"])
    code["num_coefficients"] = ret(ir["num_coefficients"])

    code["create_coordinate_finite_element"] = _create_coordinate_finite_element(ir)
    code["create_coordinate_dofmap"] = _create_coordinate_dofmap(ir)
    code["create_coordinate_mapping"] = _create_coordinate_mapping(ir)

    code["create_finite_element"] = _create_finite_element(ir)
    code["create_dofmap"] = _create_dofmap(ir)

    code["additional_includes_set"] = _additional_includes_form(ir)

    for integral_type in ufc_integral_types:
        code["max_%s_subdomain_id" % integral_type] = \
            ret(ir["max_%s_subdomain_id" % integral_type])
        code["has_%s_integrals" % integral_type] = \
            _has_foo_integrals(ir, integral_type)
        code["create_%s_integral" % integral_type] = \
            _create_foo_integral(ir, integral_type)
        code["create_default_%s_integral" % integral_type] = \
            _create_default_foo_integral(ir, integral_type)

    # Postprocess code
    _postprocess_code(code, parameters)

    return code

#--- Code generation for non-trivial functions ---


def _value_dimension(ir):
    "Generate code for value_dimension."
    ret = format["return"]
    axis = format["argument axis"]
    f_int = format["int"]

    if ir == ():
        return ret(1)
    return format["switch"](axis, [ret(n) for n in ir], ret(f_int(0)))


def _needs_mesh_entities(ir):
    """
    Generate code for needs_mesh_entities. ir is a list of num dofs
    per entity.
    """
    ret = format["return"]
    boolean = format["bool"]
    dimension = format["argument dimension"]

    return format["switch"](dimension, [ret(boolean(c)) for c in ir],
                            ret(boolean(False)))


def _global_dimension(ir):
    """Generate code for global_dimension. ir[0] is a list of num dofs per
    entity."""

    num_dofs = ir[0]
    component = format["component"]
    entities = format["num entities"]
    dimension = format["inner product"]([format["int"](d) for d in num_dofs],
                                        [component(entities, d)
                                         for d in range(len(num_dofs))])
    # Handle global "elements" if any
    if ir[1]:
        dimension = format["add"]([dimension, format["int"](ir[1])])
        try:
            dimension = format["int"](eval(dimension))
        except:
            pass

    code = format["return"](dimension)

    return code


def _tabulate_facet_dofs(ir):
    "Generate code for tabulate_facet_dofs."

    assign = format["assign"]
    component = format["component"]
    dofs = format["argument dofs"]
    cases = ["\n".join(assign(component(dofs, i), dof)
                       for (i, dof) in enumerate(facet))
             for facet in ir]
    return format["switch"](format["facet"](None), cases)


def _tabulate_dofs(ir):
    "Generate code for tabulate_dofs."

    # Prefetch formats
    add = format["addition"]
    iadd = format["iadd"]
    multiply = format["multiply"]
    assign = format["assign"]
    component = format["component"]
    entity_index = format["entity index"]
    num_entities_format = format["num entities"]
    unsigned_int = format["uint declaration"]
    dofs_variable = format["argument dofs"]

    if ir is None:
        return assign(component(dofs_variable, 0), 0)

    # Extract representation
    (dofs_per_element, num_dofs_per_element, need_offset, fakes) = ir

    # Declare offset if needed
    code = []
    offset_name = "0"
    if need_offset:
        offset_name = "offset"
        code.append(format["declaration"](unsigned_int, offset_name, 0))

    # Generate code for each element
    i = 0
    for (no, num_dofs) in enumerate(dofs_per_element):

        # Handle fakes (Space of reals)
        if fakes[no] and num_dofs_per_element[no] == 1:
            code.append(assign(component(dofs_variable, i), offset_name))
            if offset_name != "0":
                code.append(iadd(offset_name, 1))
            i += 1
            continue

        # Generate code for each degree of freedom for each dimension
        for (dim, num) in enumerate(num_dofs):

            # Ignore if no dofs for this dimension
            if not num[0]:
                continue

            for (k, dofs) in enumerate(num):
                v = multiply([len(num[k]), component(entity_index, (dim, k))])
                for (j, dof) in enumerate(dofs):
                    value = add([offset_name, v, j])
                    code.append(assign(component(dofs_variable, dof + i), value))

            # Update offset corresponding to mesh entity:
            if need_offset:
                addition = multiply([len(num[0]),
                                     component(num_entities_format, dim)])
                code.append(iadd("offset", addition))

        i += num_dofs_per_element[no]

    return "\n".join(code)


def _tabulate_dof_coordinates(ir):
    "Generate code for tabulate_dof_coordinates."

    # Raise error if tabulate_dof_coordinates is ill-defined
    if not ir:
        msg = "tabulate_dof_coordinates is not defined for this element"
        return format["exception"](msg)

    # Extract formats:
    inner_product = format["inner product"]
    component = format["component"]
    precision = format["float"]
    assign = format["assign"]
    f_x = format["coordinate_dofs"]
    coordinates = format["argument coordinates"]

    # Extract coordinates and cell dimension
    gdim = ir["gdim"]
    tdim = ir["tdim"]

    # Aid mapping points from reference to physical element
    coefficients = affine_weights(tdim)

    # Generate code for each point and each component
    code = []
    for (i, coordinate) in enumerate(ir["points"]):

        w = coefficients(coordinate)
        for j in range(gdim):
            # Compute physical coordinate
            coords = [component(f_x(), (k * gdim + j,)) for k in range(tdim + 1)]
            value = inner_product(w, coords)

            # Assign coordinate
            code.append(assign(component(coordinates, (i * gdim + j)), value))

    return "\n".join(code)


def _tabulate_entity_dofs(ir):
    "Generate code for tabulate_entity_dofs."

    # Extract variables from ir
    entity_dofs, num_dofs_per_entity = ir

    # Prefetch formats
    assign = format["assign"]
    component = format["component"]
    f_d = format["argument dimension"]
    f_i = format["argument entity"]
    dofs = format["argument dofs"]

    # Add check that dimension and number of mesh entities is valid
    dim = len(num_dofs_per_entity)
    excpt = format["exception"]("%s is larger than dimension (%d)"
                                % (f_d, dim - 1))
    code = [format["if"]("%s > %d" % (f_d, dim - 1), excpt)]

    # Generate cases for each dimension:
    all_cases = ["" for d in range(dim)]
    for d in range(dim):

        # Ignore if no entities for this dimension
        if num_dofs_per_entity[d] == 0:
            continue

        # Add check that given entity is valid:
        num_entities = len(entity_dofs[d].keys())
        excpt = format["exception"]("%s is larger than number of entities (%d)"
                                    % (f_i, num_entities - 1))
        check = format["if"]("%s > %d" % (f_i, num_entities - 1), excpt)

        # Generate cases for each mesh entity
        cases = ["\n".join(assign(component(dofs, j), dof)
                           for (j, dof) in enumerate(entity_dofs[d][entity]))
                 for entity in range(num_entities)]

        # Generate inner switch with preceding check
        all_cases[d] = "\n".join([check, format["switch"](f_i, cases)])

    # Generate outer switch
    code.append(format["switch"](f_d, all_cases))

    return "\n".join(code)


def _tabulate_entity_closure_dofs(ir):
    "Generate code for tabulate_entity_closure_dofs."

    # Extract variables from ir
    entity_closure_dofs, entity_dofs, num_dofs_per_entity = ir

    # Prefetch formats
    assign = format["assign"]
    component = format["component"]
    f_d = format["argument dimension"]
    f_i = format["argument entity"]
    dofs = format["argument dofs"]

    # Add check that dimension and number of mesh entities is valid
    dim = len(num_dofs_per_entity)
    excpt = format["exception"]("%s is larger than dimension (%d)"
                                % (f_d, dim - 1))
    code = [format["if"]("%s > %d" % (f_d, dim - 1), excpt)]

    # Generate cases for each dimension:
    all_cases = ["" for d in range(dim)]
    for d in range(dim):
        num_entities = len(entity_dofs[d])

        # Add check that given entity is valid:
        excpt = format["exception"]("%s is larger than number of entities (%d)"
                                    % (f_i, num_entities - 1))
        check = format["if"]("%s > %d" % (f_i, num_entities - 1), excpt)

        # Generate cases for each mesh entity
        cases = []
        for entity in range(num_entities):
            assignments = [assign(component(dofs, j), dof)
                           for (j, dof) in enumerate(entity_closure_dofs[(d, entity)])]
            cases.append("\n".join(assignments))

        # Generate inner switch with preceding check
        all_cases[d] = "\n".join([check, format["switch"](f_i, cases)])

    # Generate outer switch
    code.append(format["switch"](f_d, all_cases))

    return "\n".join(code)


#--- Utility functions ---

def _create_switch(arg, classnames, factory=False):
    "Generate code for create_<bar>(arg) returning new <classnames[arg]>."
    ret = format["return"]
    if factory:
        create = format["create factory"]
    else:
        create = format["create foo"]
    numbers = list(range(len(classnames)))
    cases = [ret(create(name)) for name in classnames]
    default = ret(0)
    return format["switch"](arg, cases, default=default, numbers=numbers)


def _create_coordinate_finite_element(ir):
    ret = format["return"]
    if ir["jit"]:
        create = format["create factory"]
    else:
        create = format["create foo"]
    classnames = ir["create_coordinate_finite_element"]
    assert len(classnames) == 1  # list of length 1 until we support multiple domains
    return ret(create(classnames[0]))


def _create_coordinate_dofmap(ir):
    ret = format["return"]
    if ir["jit"]:
        create = format["create factory"]
    else:
        create = format["create foo"]
    classnames = ir["create_coordinate_dofmap"]
    assert len(classnames) == 1  # list of length 1 until we support multiple domains
    return ret(create(classnames[0]))


def _create_coordinate_mapping(ir):
    ret = format["return"]
    if ir["jit"]:
        create = format["create factory"]
    else:
        create = format["create foo"]
    classnames = ir["create_coordinate_mapping"]
    assert len(classnames) == 1  # list of length 1 until we support multiple domains
    # return ret(create(classnames[0]))
    return ret("nullptr")  # FIXME: Disabled until we generate a functional class (work in progress)


def _create_finite_element(ir):
    f_i = format["argument sub"]
    classnames = ir["create_finite_element"]
    return _create_switch(f_i, classnames, ir["jit"])


def _create_dofmap(ir):
    f_i = format["argument sub"]
    classnames = ir["create_dofmap"]
    return _create_switch(f_i, classnames, ir["jit"])


def _create_sub_element(ir):
    f_i = format["argument sub"]
    classnames = ir["create_sub_element"]
    return _create_switch(f_i, classnames, ir["jit"])


def _create_sub_dofmap(ir):
    f_i = format["argument sub"]
    classnames = ir["create_sub_dofmap"]
    return _create_switch(f_i, classnames, ir["jit"])


def _has_foo_integrals(ir, integral_type):
    ret = format["return"]
    b = format["bool"]
    i = ir["has_%s_integrals" % integral_type]
    return ret(b(i))


def _create_foo_integral(ir, integral_type):
    "Generate code for create_<foo>_integral."
    ret = format["return"]
    create = format["create foo"]
    f_i = format["argument subdomain"]
    subdomain_ids, classnames = ir["create_%s_integral" % integral_type]
    cases = [ret(create(name)) for name in classnames]
    default = ret(0)
    return format["switch"](f_i, cases, default=default, numbers=subdomain_ids)


def _create_default_foo_integral(ir, integral_type):
    "Generate code for create_default_<foo>_integral."
    ret = format["return"]
    classname = ir["create_default_%s_integral" % integral_type]
    if classname is None:
        return ret(0)
    else:
        create = format["create foo"]
        return ret(create(classname))


def _postprocess_code(code, parameters):
    "Postprocess generated code."
    _indent_code(code)
    _remove_code(code, parameters)


def _indent_code(code):
    "Indent code that should be indented."
    for key in code:
        if key not in ("classname", "members", "constructor_arguments",
                       "initializer_list", "additional_includes_set",
                       "class_type"):
            code[key] = indent(code[key], 4)


def _remove_code(code, parameters):
    "Remove code that should not be generated."
    for key in code:
        flag = "no-" + key
        if flag in parameters and parameters[flag]:
            msg = "// Function %s not generated (compiled with -f%s)" \
                  % (key, flag)
            code[key] = format["exception"](msg)
