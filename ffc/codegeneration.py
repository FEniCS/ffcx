"""
Compiler stage 4: Code generation
---------------------------------

This module implements the generation of C++ code for the body of each
UFC function from an (optimized) intermediate representation (OIR).
"""

__author__ = "Anders Logg (logg@simula.no) and friends"
__date__ = "2009-12-16"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-22

# FFC modules
from ffc.log import begin, end, debug_code
from ffc.cpp import format, indent

# FFC code generation modules
from ffc.evaluatebasis import _evaluate_basis

# FFC specialized code generation modules
#from ffc.quadrature import generate_quadrature_integrals
from ffc.tensor import generate_tensor_integrals

def generate_code(ir, options):
    "Generate code from intermediate representation."

    begin("Compiler stage 4: Generating code")

    # Extract representations
    ir_form, ir_elements, ir_dofmaps, ir_integrals = ir

    # Generate code for elements, dofmaps, forms and integrals
    code_elements  = [generate_element_code(ir, options) for ir in ir_elements]
    code_dofmaps   = [generate_dofmap_code(ir, options) for ir in ir_dofmaps]
    code_integrals = generate_integrals_code(ir_integrals, options)
    code_form      = generate_form_code(ir_form, options)

    end()

    return code_form, code_elements, code_dofmaps, code_integrals

    # Generate common code like finite elements, dof map etc.
    #common_code = generate_common_code(form_data, format)

    # Generate code for integrals
    #codes = []
    #for (i, CodeGenerator) in enumerate(CodeGenerators):
    #    code_generator = CodeGenerator(options)
    #    codes.append(code_generator.generate_integrals(representations[i], format))

    # Loop all subdomains of integral types and combine code
    #combined_code = generate_combined_code(codes, form_data, prefix, format)

    # Collect generated code
    code = {}
    #code.update(common_code)
    #code.update(combined_code)

    end()
    return code

def generate_element_code(ir, options):
    "Generate code for finite element from intermediate representation."

    # Prefetch formatting to speedup code generation
    ret = format["return"]

    # Generate code
    code = {}
    code["classname"] = "FooFiniteElement"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""
    code["signature"] = ret('"%s"' % ir["signature"])
    code["cell_shape"] = ""
    code["space_dimension"] = ret(ir["space_dimension"])
    code["value_rank"] = ret(ir["value_rank"])
    code["value_dimension"] = ""
    code["evaluate_basis"] = _evaluate_basis(ir["evaluate_basis"])
#    print "CODE:\n", code["evaluate_basis"]
    code["evaluate_basis_all"] = ""
    code["evaluate_basis_derivatives"] = ""
    code["evaluate_basis_derivatives_all"] = ""
    code["evaluate_dof"] = ""
    code["evaluate_dofs"] = ""
    code["interpolate_vertex_values"] = ""
    code["num_sub_elements"] = ret(ir["num_sub_elements"])
    code["create_sub_element"] = ""

    # Postprocess code
    _postprocess_code(code, options)

    debug_code(code, "finite_element")

    return code

def generate_dofmap_code(ir, options):
    "Generate code for dofmap from intermediate representation."

    # Prefetch formatting to speedup code generation
    ret = format["return"]
    switch = format["switch"]

    # Generate code
    code = {}
    code["classname"] = "FooDofMap"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""
    code["signature"] = ret(ir["signature"])
    code["needs_mesh_entities"] = _needs_mesh_entities(ir["needs_mesh_entities"])
    code["init_mesh"] = _init_mesh(ir["init_mesh"])
    code["init_cell"] = "// Do nothing"
    code["init_cell_finalize"] = "// Do nothing"
    code["global_dimension"] = ret("__global_dimension")
    code["local_dimension"] = ret(ir["local_dimension"])
    code["max_local_dimension"] = ret(ir["max_local_dimension"])
    code["geometric_dimension"] = ret(ir["geometric_dimension"])
    code["num_facet_dofs"] = ret(ir["num_facet_dofs"])
    code["num_entity_dofs"] = "// Marie does not know what this should return // AL: Should return number of dofs associated with the entity of a given topological dimension, so 0 --> 1, 1 --> 2, 2 --> 1 for P3 triangles"
    code["tabulate_dofs"] = ""
    code["tabulate_facet_dofs"] = ""
    code["tabulate_entity_dofs"] = ""
    code["tabulate_coordinates"] = ""
    code["num_sub_dof_maps"] = ret(ir["num_sub_dof_maps"])
    code["create_sub_dof_map"] = ""

    # Postprocess code
    _postprocess_code(code, options)

    debug_code(code, "dofmap")

    return code

def generate_integrals_code(ir, options):
    "Generate code for integrals from intermediate representation."

    # FIXME: Handle multiple representations here

    code = generate_tensor_integrals(ir, options)

    print code

    return code

def generate_form_code(ir, options):
    "Generate code for form from intermediate representation."

    # Prefetch formatting to speedup code generation
    ret = format["return"]

    # Generate code
    code = {}
    code["classname"] = "FooForm"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""
    code["signature"] = ret('"%s"' % ir["signature"])
    code["rank"] = ret(ir["rank"])
    code["num_coefficients"] = ret(ir["num_coefficients"])
    code["num_cell_integrals"] = ret(ir["num_cell_integrals"])
    code["num_exterior_facet_integrals"] = ret(ir["num_exterior_facet_integrals"])
    code["num_interior_facet_integrals"] = ret(ir["num_interior_facet_integrals"])
    code["create_finite_element"] = ""
    code["create_dof_map"] = ""
    code["create_cell_integral"] = ""
    code["create_exterior_facet_integral"] = ""
    code["create_interior_facet_integral"] = ""

    # Postprocess code
    _postprocess_code(code, options)

    debug_code(code, "form")

    return code

def _needs_mesh_entities(num_dofs_per_entity):
    "Generate code for needs_mesh_entities."
    return format["switch"]("d", [format["return"](format["bool"](c))
                                  for c in num_dofs_per_entity])

def _init_mesh(num_dofs_per_entity):
    "Generate code for init_mesh."
    terms = [format["multiply"](["%d" % num, "m.num_entities[%d]" % dim])
             for (dim, num) in enumerate(num_dofs_per_entity)]
    dimension = format["add"](terms)
    return "__global_dimension = %s;\n return false;" % dimension

def _postprocess_code(code, options):
    "Postprocess generated code."
    _indent_code(code)
    _remove_code(code, options)

def _indent_code(code):
    "Indent code that should be indented."
    for key in code:
        if not key in ("classname", "members"):
            code[key] = indent(code[key], 4)

def _remove_code(code, options):
    "Remove code that should not be generated."
    for key in code:
        flag = "no-" + key
        if flag in options and options[flag]:
            msg = "// Function %s not generated (compiled with -f%s)" % (key, flag)
            code[key] = format["exception"](msg)
