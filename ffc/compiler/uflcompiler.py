"""This is the compiler, acting as the main interface for compilation
of forms and breaking the compilation into several sequential stages:

   0. language        -  expressing the form in the form language (UFL)
   1. analysis        -  simplifying and preprocessing the form (UFL)
   2. representation  -  computing a representation of the form (FIAT/FFC)
   3. optimization    -  optimizing the form representation (FErari)
   4. codegeneration  -  generating code according to a format (FFC)
   5. format          -  writing the generated code to file (FFC -> UFC)

"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-05 -- 2009-01-13"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.
# Modified by Dag Lindbo, 2008.

# UFL modules
from ufl.classes import Form, FiniteElementBase, Measure, Integral
from ufl.algorithms import FormData, is_multilinear, validate_form, extract_monomials, extract_quadrature_order, estimate_quadrature_order


# FFC common modules
from ffc.common.log import debug, info, warning, error, begin, end, set_level, INFO
from ffc.common.constants import FFC_OPTIONS

# FFC fem modules
#from ffc.fem.finiteelement import *
#from ffc.fem.vectorelement import *
#from ffc.fem.mixedelement import *
#from ffc.fem.mixedfunctions import *

# FIXME: Remove this
from analysis.formdata import FormData as FFCFormData

# FFC form representation modules
from representation.tensor import UFLTensorRepresentation
from representation.quadrature import UFLQuadratureRepresentation

# FFC code generation modules
#from codegeneration.tensor import *
#from codegeneration.quadrature import *
from codegeneration.common.codegenerator import CodeGenerator
from codegeneration.quadrature import UFLQuadratureGenerator
#from codegeneration.common.finiteelement import *
#from codegeneration.common.dofmap import *

# FFC format modules
from format import uflufcformat
#from format import dolfinformat

def compile(objects, prefix="Form", options=FFC_OPTIONS):
    "Compile the given forms and/or elements."

    # Set log level
    set_level(INFO)

    # Check options
    _check_options(options)

    # Extract forms and elements
    (forms, elements) = _extract_objects(objects)
    if len(forms) + len(elements) == 0:
        info("No forms or elements specified, nothing to do.")
        return

    # Compile forms
    #(form_data, form_representation) = _compile_forms(forms, prefix, options)
    if len(forms) > 0:
        compile_forms(forms, prefix, options)

    # Compile elements
    # FIXME/TODO: Do we still need this?
    #if len(elements) > 0:
    #    compile_elements(elements, prefix, options)

    #    return (form_data, form_representation)

def compile_forms(forms, prefix, options):
    "Compile the given forms."
    # FIXME: Make format a class

    # Choose format
    format = uflufcformat
    format.init(options)

    # Iterate over forms for stages 1 - 4
    generated_forms = []
#    form_datas = []
#    form_representations = []
    for form in forms:

        # Handle all metadata of integrals
        # TODO: Improve algorithm in UFL that returns the quad_order of
        # an integral.
        # NOTE: As a result of handling the metadata, some measures
        # might become equal which means that integrals can be grouped.
        # The handle_metadatas() must therefore be called before analyze_form()
        # like it is now.
        form = handle_metadatas(form, options)

        # Compiler stage 1: analyze form
        form_data = analyze_form(form)
        #form_datas += [form_data]

        # Compiler stage 2: compute form representation
        tensor_representation, quadrature_representation =\
            compute_form_representation(form_data, options)

        # Compiler stage 3: optimize form representation
        # TODO: Switch this back on? I guess the argument should only be the
        # tensor_representation since it only applies to this.
        #optimize_form_representation(form)

        # Compiler stage 4: generate form code
        #form_code = generate_form_code(form_data, form_representation, options["representation"], format.format)
        ffc_form_data = FFCFormData(None, ufl_form_data=form_data)
        form_code = generate_form_code(ffc_form_data, tensor_representation,\
                          quadrature_representation, format.format)

        # Add to list of codes
        generated_forms += [(form_code, ffc_form_data)]

    # Compiler stage 5: format code
    format_code(generated_forms, prefix, format, options)

    return
#    return (form_datas, form_representations)

def compile_elements(elements, prefix="Element", options=FFC_OPTIONS):
    "Compile the given elements for the given language"

    # Check input
    if len(forms) == 0:
        return

    # Compiler stage 1: analyze form
    begin("Compiler stage 1: Analyzing elements")
    element_data = ElementData(elements)
    end()

    # Go directly to stage 4, code generation
    begin("Compiler stage 2-3: Nothing to do")
    end()
    begin("Compiler stage 4: Generating code")

    # Choose format
    format = __choose_format(options)
    format.init(options)

    # Choose code generator
    CodeGenerator = __choose_code_generator(options["representation"])
    code_generator = CodeGenerator()

    # Generate code
    element_code = code_generator.generate_element_code(element_data, format.format)

    # Compiler stage 5: format code
    end()
    begin("Compiler stage 5: Formatting code")
    format.write([(element_code, element_data)], prefix, options)
    end()

def handle_metadatas(form, options):
    "Handle metadata of all integrals"

    # TODO: Is this the best way of doing this?
    # Create new and empty Form
    form_new = Form([])

    # Loop all integrals and create new forms. Add these forms such that
    # integrals which might have become equal (due to metadata handling)
    # are grouped.
    for integral in form.cell_integrals():
        form_new += Form([handle_metadata(integral, options)])
    for integral in form.exterior_facet_integrals():
        form_new += Form([handle_metadata(integral, options)])
    for integral in form.interior_facet_integrals():
        form_new += Form([handle_metadata(integral, options)])
#    print "new form: ", form_new

    return form_new

def handle_metadata(integral, options):
    "Handle metadata of one integral"

    # Get the old metadata
    metadata_old = integral.measure().metadata()

    # Set default values for representation and quadrature_order
    representation = options["representation"]
    quadrature_order = options["quadrature_order"]

    # If we have a metadata loop options
    if metadata_old:
        for k, v in metadata_old.items():
            if k == "ffc":
                if "representation" in v:
                    # Get representation and check that it is valid
                    representation = v.pop("representation")
                    if not representation in ["tensor", "quadrature", "automatic"]:
                        error("Unrecognized representation '%s', must be one of: ['tensor', 'quadrature', 'automatic']" % representation)
                    # If we still have some options, display a warning
                    if v:
                        warning("Following options are not supported: " + ", ".join([str(o) for o in v]))
            elif k == "quadrature_order":
                quadrature_order = v
            else:
                warning("Unrecognized option %s for metadata" % k)

    # Automatically select representation based on operation estimate
    if representation == "automatic":
        representation = auto_select_representation(integral)

    # TODO: If quadrature_order can't be converted to an int it will fail
    if quadrature_order != None:
        quadrature_order = int(quadrature_order)
    else:
        # TODO: Improve algorithm in UFL
        quadrature_order = max(extract_quadrature_order(integral),\
                               estimate_quadrature_order(integral))

    # Create the new consistent metadata and Measure
    metadata_new = {"quadrature_order":quadrature_order, "ffc":{"representation":representation}}
    measure_new = integral.measure().reconstruct(metadata=metadata_new)

    # Create and return new Integral
    return Integral(integral.integrand(), measure_new)

def auto_select_representation(integral):
    "Automatically select the best representation"

    # FIXME: Implement this
    return "quadrature"

def analyze_form(form):
    "Compiler stage 1: analyze form"
    begin("Compiler stage 1: Analyzing form")
    validate_form(form)
    form_data = FormData(form)
    info(str(form_data))
    end()
    return form_data

def compute_form_representation(form_data, options):
    "Compiler stage 2: Compute form representation"
    begin("Compiler stage 2: Computing form representation")

    # Choose representation
#    Representation = _choose_representation(form_data.form, options)

    # Testing monomial extraction
#    print form_data.form
#    monomial = extract_monomials(form_data.form)
#    print monomial

#    import sys
#    sys.exit(0)

    # Compute form representation for both representations
    # TODO: Switch on tensor
#    tensor = UFLTensorRepresentation(form_data)
    quadrature = UFLQuadratureRepresentation(form_data)

    end()
    return (quadrature, quadrature)
#    return (tensor, quadrature)
    
def optimize_form_representation(form):
    "Compiler stage 3: Compute optimization"
    begin("Compiler stage 3: Computing optimization")
    info("Optimization currently broken (to be fixed).")
    end()

def generate_form_code(form_data, tensor_representation, quadrature_representation, format):
    "Compiler stage 4: Generate code"
    begin("Compiler stage 4: Generating code")

    # Generate common code like finite elements, dof map etc.
    common_generator = CodeGenerator()
    code = common_generator.generate_form_code(form_data, None, format, ufl_code=True)

    # We need both genrators
#    tensor_generator = UFLTensorGenerator()
    tensor_generator = UFLQuadratureGenerator()
    quadrature_generator = UFLQuadratureGenerator()

    # Generate code for integrals using quadrature
    quadrature_code = quadrature_generator.generate_cell_integrals(quadrature_representation, format)
    quadrature_code.update(quadrature_generator.generate_exterior_facet_integrals(quadrature_representation, format))
    quadrature_code.update(quadrature_generator.generate_interior_facet_integrals(quadrature_representation, format))

    # Mock tensor code
    tensor_code = {}
#    tensor_code.update(tensor_generator.generate_cell_integrals(quadrature_representation, format))
#    tensor_code.update(tensor_generator.generate_exterior_facet_integrals(quadrature_representation, format))
#    tensor_code.update(tensor_generator.generate_interior_facet_integrals(quadrature_representation, format))

#    code["quadrature"].update(quadrature_code)
#    code.update(quadrature_code)

    ## A first try at combining code from the two generators
    ## Could probably be moved somewhere else

    # Get any kind of code for reseting the element tensor, it just needs to be
    # generated once by the codegenerators
    reset_code = (quadrature_generator.reset_code or tensor_generator.reset_code)
    reset_code_restricted = (quadrature_generator.reset_code_restricted or tensor_generator.reset_code_restricted)

    # Loop all subdomains of integral types and combine code
    # TODO: Only cell_integrals are handled in uflufcformat.py and split_implementation
    # is not handled
    for i in range(form_data.num_cell_integrals):
        combine_code(code, tensor_code, quadrature_code, ("cell_integral", i), reset_code, False)
    for i in range(form_data.num_exterior_facet_integrals):
        combine_code(code, tensor_code, quadrature_code, ("exterior_facet_integral", i), reset_code, True)
    for i in range(form_data.num_interior_facet_integrals):
        combine_code(code, tensor_code, quadrature_code, ("interior_facet_integral", i), reset_code_restricted, True)

    end()
    return code

def combine_code(code, tensor_code, quadrature_code, key, reset_code, facet_integral):
    "Combine the code from the two code generators"

    # If subdomain has both representations then combine them
    if key in tensor_code and key in quadrature_code:
        code[key] = {("tabulate_tensor_tensor"):tensor_code[key],
                     ("tabulate_tensor_quadrature"):quadrature_code[key],
                     "reset_tensor": reset_code}
    # Handle code from tensor generator
    elif key in tensor_code:
        # Add reset code to tabulate_tensor code
        val = tensor_code[key]["tabulate_tensor"]
        # Check if we have a tuple (common, cases) for facet integrals
        if facet_integral:
            val = (reset_code + val[0], val[1])
        else:
            val = reset_code + val
        tensor_code[key]["tabulate_tensor"] = val
        code[key] = tensor_code[key]

    # Handle code from quadrature generator
    elif key in quadrature_code:
        # Add reset code to tabulate_tensor code
        val = quadrature_code[key]["tabulate_tensor"]
        # Check if we have a tuple (common, cases) for facet integrals
        if facet_integral:
            val = (reset_code + val[0], val[1])
        else:
            val = reset_code + val
        quadrature_code[key]["tabulate_tensor"] = val
        code[key] = quadrature_code[key]

    # If we reach this level it means that no code has been generated
    # for the given subdomain so we need to add the reset code
    # NOTE: If we were sure that all assemblers would reset the local
    # tensor before calling tabulate_tensor this wouldn't be needed
    else:
        if facet_integral:
            code[key] = {"tabulate_tensor": (reset_code, []), "members": []}
        else:
            code[key] = {"tabulate_tensor": reset_code, "members": []}


def format_code(generated_forms, prefix, format, options):
    "Compiler stage 5: Format code"
    begin("Compiler stage 5: Formatting code")
    format.write(generated_forms, prefix, options)
    end()

def _check_options(options):
    "Check that options are valid"
    if options["optimize"]:
        warning("Optimization unavailable (will return in a future version).")
    if options["blas"]:
        warning("BLAS mode unavailable (will return in a future version).")
    if options["quadrature_points"]:
        warning("The option 'quadrature_points' is only available for standard FFC forms (not UFL forms), will be ignored.")

def _extract_objects(objects):
    "Extract forms and elements from list of objects."

    # Check that we get a list
    if not isinstance(objects, list):
        objects = [objects]

    # Check each object
    forms = []
    elements = []
    for object in objects:
        if isinstance(object, Form):
            forms.append(object)
        elif isinstance(object, FiniteElementBase):
            elements.append(object)
        elif not object is None:
            error("Not a form: " + str(object))

    # Only compile element(s) when there are no forms
    if len(forms) > 0 and len(elements) > 0:
        elements = []

    return (forms, elements)

def _choose_format(options):
    "Choose output format."

    # FIXME: Make format a class (since we call init())

    language = options["language"]
    if language.lower() == "ufc":
        return ufcformat
    elif language.lower() == "dolfin":
        return dolfinformat
    else:
        raise RuntimeError, "Don't know how to compile code for language \"%s\"." % language

def _choose_representation(form, options):
    "Choose form representation"

    option = options["representation"]
    if option == "tensor":

        # Check if form is multilinear
        info("Checking if form is multilinear...")
        if is_multilinear(form):
            info("yes\n")
            return UFLTensorRepresentation
        else:
            info("no\n")
            warning("Form is is not multilinear, using quadrature representation")
            return UFLQuadratureRepresentation
        
    elif option == "quadrature":
        return UFLQuadratureRepresentation

    else:
        raise RuntimeError, 'Unknown form representation: "%s"' % option

def __choose_code_generator(form_representation):
    "Choose code generator"

    if form_representation == "tensor":
        return TensorGenerator
    else:
        return UFLQuadratureGenerator
    "Choose code generator"

    if form_representation == "tensor":
        return TensorGenerator
    else:
        return UFLQuadratureGenerator
