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

# Last changed: 2010-01-07

# FFC modules
from ffc.log import begin, end, debug_code
from ffc.cpp import format, indent

# FFC code generation modules
from ffc.evaluatebasis import _evaluate_basis
from ffc.evaluatedof import _evaluate_dof, _evaluate_dofs, affine_weights
from ffc.codesnippets import jacobian

# FFC specialized code generation modules
#from ffc.quadrature import generate_quadrature_integrals
from ffc.tensor import generate_tensor_integrals

# FIXME: Temporary
not_implemented = "// Not implemented, please fix me!"

def generate_code(ir, prefix, options):
    "Generate code from intermediate representation."

    begin("Compiler stage 4: Generating code")

    # Extract representations
    ir_form, ir_elements, ir_dofmaps, ir_integrals = ir

    # Generate code for elements
    code_elements = [generate_element_code(i, ir, prefix, options)
                     for (i, ir) in enumerate(ir_elements)]

    # Geneate code for dofmaps
    code_dofmaps = [generate_dofmap_code(i, ir, prefix, options)
                    for (i, ir) in enumerate(ir_dofmaps)]

    # Generate code for integrals
    code_integrals = generate_integrals_code(ir_integrals, prefix, options)

    # Generate code for form
    code_form = generate_form_code(ir_form, prefix, options)

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

def generate_element_code(i, ir, prefix, options):
    "Generate code for finite element from intermediate representation."

    # Prefetch formatting to speedup code generation
    ret = format["return"]
    do_nothing = format["do nothing"]

    # Generate code
    code = {}
    code["classname"] = prefix.lower() + "_finite_element_" + str(i)
    code["members"] = ""
    code["constructor"] = do_nothing
    code["destructor"] = do_nothing
    code["signature"] = ret('"%s"' % ir["signature"])
    code["cell_shape"] = ret("ufc:%s" % ir["cell_shape"])
    code["space_dimension"] = ret(ir["space_dimension"])
    code["value_rank"] = ret(ir["value_rank"])
    code["value_dimension"] = _value_dimension(ir["value_dimension"])
    code["evaluate_basis"] = _evaluate_basis(ir["evaluate_basis"])
    code["evaluate_basis_all"] = not_implemented
    code["evaluate_basis_derivatives"] = not_implemented
    code["evaluate_basis_derivatives_all"] = not_implemented
    code["evaluate_dof"] = _evaluate_dof(ir["evaluate_dof"])
    code["evaluate_dofs"] = _evaluate_dofs(ir["evaluate_dofs"])
    code["interpolate_vertex_values"] = _interpolate_vertex_values(ir["interpolate_vertex_values"])
    code["num_sub_elements"] = ret(ir["num_sub_elements"])
    code["create_sub_element"] = _create_foo(ir["create_sub_element"], prefix, "finite_element")

    # Postprocess code
    _postprocess_code(code, options)

    debug_code(code, "finite_element")

    return code

def generate_dofmap_code(i, ir, prefix, options):
    "Generate code for dofmap from intermediate representation."

    # Prefetch formatting to speedup code generation
    ret = format["return"]
    do_nothing = format["do nothing"]

    # Generate code
    code = {}
    code["classname"] = prefix.lower() + "_dof_map_" + str(i)
    code["members"] = ""
    code["constructor"] = do_nothing
    code["destructor"] = do_nothing
    code["signature"] = ret('"%s"' % ir["signature"])
    code["needs_mesh_entities"] = _needs_mesh_entities(ir["needs_mesh_entities"])
    code["init_mesh"] = _init_mesh(ir["init_mesh"])
    code["init_cell"] = do_nothing
    code["init_cell_finalize"] = do_nothing
    code["global_dimension"] = ret("__global_dimension")
    code["local_dimension"] = ret(ir["local_dimension"])
    code["max_local_dimension"] = ret(ir["max_local_dimension"])
    code["geometric_dimension"] = ret(ir["geometric_dimension"])
    code["num_facet_dofs"] = ret(ir["num_facet_dofs"])
    code["num_entity_dofs"] = format["switch"]("d", [ret(num) for num in ir["num_entity_dofs"]])
    code["tabulate_dofs"] = _tabulate_dofs(ir["tabulate_dofs"])
    code["tabulate_facet_dofs"] = _tabulate_facet_dofs(ir["tabulate_facet_dofs"])
    code["tabulate_entity_dofs"] = not_implemented # Marie doesn't know what this function should do
    code["tabulate_coordinates"] = _tabulate_coordinates(ir["tabulate_coordinates"])
    code["num_sub_dof_maps"] = ret(ir["num_sub_dof_maps"])
    code["create_sub_dof_map"] = _create_foo(ir["create_sub_dof_map"], prefix, "dof_map")

    # Postprocess code
    _postprocess_code(code, options)

    debug_code(code, "dofmap")

    return code

def generate_integrals_code(ir, prefix, options):
    "Generate code for integrals from intermediate representation."

    # FIXME: Handle multiple representations here

    code = generate_tensor_integrals(ir, options)

    # Postprocess code
    for integral_type_code in code:
        for sub_domain_code in integral_type_code:
            _postprocess_code(sub_domain_code, options)

    return code

def generate_form_code(ir, prefix, options):
    "Generate code for form from intermediate representation."

    # Prefetch formatting to speedup code generation
    ret = format["return"]
    do_nothing = format["do nothing"]

    # Generate code
    code = {}
    code["classname"] = prefix.lower() + "_form"
    code["members"] = ""
    code["constructor"] = do_nothing
    code["destructor"] = do_nothing
    code["signature"] = ret('"%s"' % ir["signature"])
    code["rank"] = ret(ir["rank"])
    code["num_coefficients"] = ret(ir["num_coefficients"])
    code["num_cell_integrals"] = ret(ir["num_cell_integrals"])
    code["num_exterior_facet_integrals"] = ret(ir["num_exterior_facet_integrals"])
    code["num_interior_facet_integrals"] = ret(ir["num_interior_facet_integrals"])
    code["create_finite_element"] = _create_foo(ir["create_finite_element"], prefix, "finite_element")
    code["create_dof_map"] = _create_foo(ir["create_dof_map"], prefix, "dof_map")
    code["create_cell_integral"] = _create_foo(ir["create_cell_integral"], prefix, "cell_integral")
    code["create_exterior_facet_integral"] = _create_foo(ir["create_exterior_facet_integral"], prefix, "exterior_facet_integral")
    code["create_interior_facet_integral"] = _create_foo(ir["create_interior_facet_integral"], prefix, "interior_facet_integral")

    # Postprocess code
    _postprocess_code(code, options)

    debug_code(code, "form")

    return code

#--- Code generation for non-trivial functions ---

def _value_dimension(ir):
    "Generate code for value_dimension."
    if ir == ():
    # FIXME: KBO: Use format instead of "1" and str(n)
        return format["return"]("1")
    return format["switch"]("i", [format["return"](str(n)) for n in ir])

def _needs_mesh_entities(ir):
    "Generate code for needs_mesh_entities."
    num_dofs_per_entity = ir
    return format["switch"]("d", [format["return"](format["bool"](c))
                                  for c in num_dofs_per_entity])

def _init_mesh(ir):
    "Generate code for init_mesh."
    num_dofs_per_entity = ir
    terms = [format["multiply"](["%d" % num, "m.num_entities[%d]" % dim])
             for (dim, num) in enumerate(num_dofs_per_entity)]
    dimension = format["add"](terms)
    return "__global_dimension = %s;\n return false;" % dimension

def _tabulate_facet_dofs(ir):
    "Generate code for tabulate_facet_dofs."
    return "\n".join([format["switch"]("facet", ["\n".join(["dofs[%d] = %d;" % (i, dof)
                                                            for (i, dof) in enumerate(ir[facet])])
                                                 for facet in range(len(ir))])])

def _tabulate_dofs(ir):
    "Generate code for for tabulate_dofs."

    # Note:  Not quite C++ independent, and not optimized.

    # Prefetch add and multiply
    add = format["add"]
    multiply = format["multiply"]

    # Declare offset
    code = ["unsigned int offset = 0;"]

    # Total dof counter
    i = 0

    # Generate code for each element
    for element_ir in ir:

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

def _tabulate_coordinates(ir):

    # Raise error if tabulate_coordinates is ill-defined
    if ir is None:
        return "// Raise error here"

    # Extract formats:
    add, multiply = format["add"], format["multiply"]
    precision = format["float"]

    # Extract coordinates and cell dimension
    cell_dim = len(ir[0])

    # Aid mapping points from reference to physical element
    coefficients = affine_weights(cell_dim)

    # Generate code for each point and each component
    code = ["const double * const * x = c.coordinates;"]
    for (i, coordinate) in enumerate(ir):
        w = coefficients(coordinate)
        for j in range(cell_dim):
            value = add([multiply([precision(w[k]), "x[%d][%d]" % (k, j)])
                         for k in range(cell_dim + 1)])
            code += ["coordinates[%d][%d] = %s;" % (i, j, value)]
    return "\n".join(code)

def _interpolate_vertex_values(ir):

    code = []

    # Add code for Jacobian if necessary
    if ir["needs_jacobian"]:
        code += [jacobian(ir["cell_dim"])]

    # Extract formats
    add, multiply = format["add"], format["multiply"]
    precision = format["float"]

    # Iterate over the elements
    for (k, data) in enumerate(ir["element_data"]):

        # Extract vertex values for all basis functions
        vertex_values = data["basis_values"]

        # Create code for each vertex
        for j in range(len(vertex_values)):

            # Map basis values according to mapping
            vertex_value = [precision(v) for v in vertex_values[j]]

            # Contract basis values and coefficients
            value = add([multiply(["dof_values[%d]" % i, v])
                         for (i, v) in enumerate(vertex_value)])

            # Construct correct vertex value label
            name = "vertex_values[%d]" % j
            code += [name + " = " +  value + ";"]

    return "\n".join(code)

#--- Utility functioins ---

def _create_foo(numbers, prefix, class_name):
    "Generate code for create_<foo>."
    class_names = ["%s_%s_%d" % (prefix, class_name, i) for i in numbers]
    cases = [format["return"]("new " + name + "()") for name in class_names]
    return format["switch"]("i", cases)

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
