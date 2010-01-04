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

# Last changed: 2010-01-04

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

    return code_elements, code_dofmaps, code_integrals, code_form

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
    code["cell_shape"] = ret("ufc:%s" % ir["cell_shape"])
    code["space_dimension"] = ret(ir["space_dimension"])
    code["value_rank"] = ret(ir["value_rank"])
    code["value_dimension"] = _value_dimension(ir["value_dimension"])
    code["evaluate_basis"] = _evaluate_basis(ir["evaluate_basis"])
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

def _value_dimension(ir):
    print "ir: ", ir
    if ir == ():
        return format["return"]("1")
    return format["switch"]("i", [n for n in ir])

def generate_dofmap_code(ir, options):
    "Generate code for dofmap from intermediate representation."

    not_implemented = "// NotImplementedYet"

    # Prefetch formatting to speedup code generation
    ret = format["return"]

    # Generate code
    code = {}
    code["classname"] = "FooDofMap"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""
    code["signature"] = ret('"%s"' % ir["signature"])
    code["needs_mesh_entities"] = _needs_mesh_entities(ir["needs_mesh_entities"])
    code["init_mesh"] = _init_mesh(ir["init_mesh"])
    code["init_cell"] = "// Do nothing"
    code["init_cell_finalize"] = "// Do nothing"
    code["global_dimension"] = ret("__global_dimension")
    code["local_dimension"] = ret(ir["local_dimension"])
    code["max_local_dimension"] = ret(ir["max_local_dimension"])
    code["geometric_dimension"] = ret(ir["geometric_dimension"])
    code["num_facet_dofs"] = ret(ir["num_facet_dofs"])
    code["num_entity_dofs"] = format["switch"]("d", [ret(num) for num in ir["num_entity_dofs"]])
    code["tabulate_dofs"] = _tabulate_dofs(ir["tabulate_dofs"])
    code["tabulate_facet_dofs"] = _tabulate_facet_dofs(ir["tabulate_facet_dofs"])
    code["tabulate_entity_dofs"] = "// Marie doesn't know what this function should do."
    code["tabulate_coordinates"] = "// Marie doesn't believe in this function."
    code["num_sub_dof_maps"] = ret(ir["num_sub_dof_maps"])
    code["create_sub_dof_map"] = not_implemented

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

def _tabulate_facet_dofs(tabulate_facet_dofs_ir):
    """
    Code generation for tabulate facet dofs.
    """

    return "\n".join([format["switch"]("facet", ["\n".join(["dofs[%d] = %d;" % (i, dof)
                                                            for (i, dof) in enumerate(tabulate_facet_dofs_ir[facet])])
                                                 for facet in range(len(tabulate_facet_dofs_ir))])])

def _tabulate_dofs(tabulate_dofs_ir):
    """Code generation for tabulate dofs.

    Not quite c++ independent, and not optimized.
    """

    # Prefetch add and multiply
    add = format["add"]
    multiply = format["multiply"]

    # Declare offset
    code = ["unsigned int offset = 0;"]

    # Total dof counter
    i = 0

    # Generate code for each element
    for element_ir in tabulate_dofs_ir:

        # Representation contains number of dofs per mesh entity and
        # number of mesh entities per geometric dimension
        dofs_per_entity = element_ir["num_dofs_per_entity"]
        entities_per_dim = element_ir["entites_per_dim"]

        # For each dimension, generate code for each dof
        for (d, num_dofs) in enumerate(dofs_per_entity):
            if num_dofs > 0:
                for k in range(entities_per_dim[d]):
                    for j in range(num_dofs):

                        name = "dofs[%d]" % i
                        value = multiply(["%d" % num_dofs, format["entity index"](d, k)])
                        value = add(["offset", value, "%d" % j])
                        code += [name + " = " + value + ";"]
                        i += 1

                # Set offset corresponding to mesh entity:
                code += [format["iadd"]("offset", multiply(["%d" % num_dofs, format["num entities"](d)]))]

    return "\n".join(code)




def _postprocess_code(code, options):
    "Postprocess generated code."
    _indent_code(code)
    _remove_code(code, options)

def _indent_code(code):
    "Indent code that should be indented."
    for key in code:
        if not key in ("classname", "members"):
            print "key: ", key
            print "code: ", code[key]
            code[key] = indent(code[key], 4)

def _remove_code(code, options):
    "Remove code that should not be generated."
    for key in code:
        flag = "no-" + key
        if flag in options and options[flag]:
            msg = "// Function %s not generated (compiled with -f%s)" % (key, flag)
            code[key] = format["exception"](msg)
