"""This is the compiler, acting as the main interface for compilation
of forms and breaking the compilation into several sequential phases:

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
from ufl.classes import Form, FiniteElementBase
from ufl.algorithms import FormData, is_multilinear
from ufl.output import set_loglevel as ufl_loglevel

# FFC common modules
from ffc.common.log import debug, info, warning, error, begin, end, set_level, INFO
from ffc.common.constants import FFC_OPTIONS

# FFC fem modules
#from ffc.fem.finiteelement import *
#from ffc.fem.vectorelement import *
#from ffc.fem.mixedelement import *
#from ffc.fem.mixedfunctions import *

# FFC language modules
#from language.algebra import *

# FFC analysis modules
#from analysis.analyze import *

# FFC form representation modules
from representation.tensor import UFLTensorRepresentation
from representation.quadrature import UFLQuadratureRepresentation

# FFC code generation modules
#from codegeneration.tensor import *
#from codegeneration.quadrature import *
#from codegeneration.common.finiteelement import *
#from codegeneration.common.dofmap import *

# FFC format modules
from format import ufcformat
from format import dolfinformat

def compile(objects, prefix="Form", options=FFC_OPTIONS):
    "Compile the given forms and/or elements."

    # Set log level
    set_level(INFO)

    # Check options
    _check_options(options)

    # Extract forms and elements
    (forms, elements) = _extract_objects(objects)
    if len(forms + elements) == 0:
        info("No forms or elements specified, nothing to do.")
        return

    # Compile forms
    #(form_data, form_representation) = _compile_forms(forms, prefix, options)
    if len(forms) > 0:
        _compile_forms(forms, prefix, options)

    # Compile elements
    #_compile_elements(elements, prefix, options)
    #_compile_elements(elements, prefix, options)
#    return (form_data, form_representation)

def _compile_forms(forms, prefix, options):
    "Compile the given forms."

    # Choose format
    format = _choose_format(options)
    format.init(options)

    # Iterate over forms for stages 1 - 4
    generated_forms = []
#    form_datas = []
#    form_representations = []
    for form in forms:

        # Compiler phase 1: analyze form
        form_data = analyze_form(form)
#        form_datas += [form_data]

        # Representations on subdomains
        # FIXME: We can support different representations on individual
        # subdomains. The representation can be specified by the user on the
        # command line (global represetation for all forms in the ufl file);
        # extracted as meta data from the integrals in the ufl file (user
        # specified for individual subdomain integrals); or the representation
        # can be selected automatically by a module yet to be implemented.
        # Key thing to note is that each subdomain uses the same representation
        # for ALL terms, otherwise we will get into all kinds of troubles with
        # with respect to the codegenerators.
        # Check if specified representation is legal.
        domain_representations = {}
        # Simple generator for now
        for integral in form_data.form.cell_integrals():
            domain_representations[(integral.domain_type(), integral.domain_id())] = options["representation"]
        for integral in form_data.form.exterior_facet_integrals():
            domain_representations[(integral.domain_type(), integral.domain_id())] = options["representation"]
        for integral in form_data.form.interior_facet_integrals():
            domain_representations[(integral.domain_type(), integral.domain_id())] = options["representation"]

        print "domain_representations:\n", domain_representations

        # Compiler phase 2: compute form representation
        tensor_representation, quadrature_representation =\
            compute_form_representation(form_data, domain_representations, options)

#        form_representation = compute_form_representation(form_data, options)
#        form_representations += [form_representation]

        # Compiler phase 3: optimize form representation
        #optimize_form_representation(form)

        # Compiler phase 4: generate form code
        #form_code = generate_form_code(form_data, form_representation, options["representation"], format.format)

        # Add to list of codes
        #generated_forms += [(form_code, form_data)]

    # Compiler phase 5: format code
    #_format_code(generated_forms, prefix, format, options)

    return
#    return (form_datas, form_representations)

def __compile_elements(elements, prefix="Element", options=FFC_OPTIONS):
    "Compile the given elements for the given language"

    # Check input
    if len(forms) == 0:
        return

    

    # Compiler phase 1: analyze form
    begin("Compiler phase 1: Analyzing elements")
    element_data = ElementData(elements)
    end()

    # Go directly to phase 4, code generation
    begin("Compiler phase 2-3: Nothing to do")
    end()
    begin("Compiler phase 4: Generating code")

    # Choose format
    format = __choose_format(options)
    format.init(options)

    # Choose code generator
    CodeGenerator = __choose_code_generator(options["representation"])
    code_generator = CodeGenerator()

    # Generate code
    element_code = code_generator.generate_element_code(element_data, format.format)

    # Compiler phase 5: format code
    end()
    begin("Compiler phase 5: Formatting code")
    format.write([(element_code, element_data)], prefix, options)
    end()
    
def analyze_form(form):
    "Compiler phase 1: analyze form"
    begin("Phase 1: Analyzing form")
    form_data = FormData(form)
    info(str(form_data))
    end()
    return form_data

def compute_form_representation(form_data, domain_representations, options):
    "Compiler phase 2: compute form representation"
    begin("Compiler phase 2: Computing form representation")

    # Choose representation
#    Representation = _choose_representation(form_data.form, options)

    # Compute form representation
    # FIXME: The representations should of course only be generated for the
    # relevant subdomains
#    tensor = UFLTensorRepresentation(form_data, int(options["quadrature_points"]))
    quadrature = UFLQuadratureRepresentation(form_data, domain_representations, int(options["quadrature_points"]))

    end()
    return (quadrature, quadrature)
#    return (tensor, quadrature)
    
def optimize_form_representation(form):
    "Compiler phase 3: optimize form representation"
    begin("Compiler phase 3: Computing optimization")

    info("Optimization currently broken (to be fixed).")

    end()

def generate_form_code(form_data, form_representation, representation, format):
    "Compiler phase 4: generate code"
    begin("Compiler phase 4: Generating code")

    # Choose code generator
    CodeGenerator = __choose_code_generator(representation)

    # Generate code
    code_generator = CodeGenerator()
    code = code_generator.generate_form_code(form_data, form_representation, format)
        
    end()
    return code

def format_code(generated_forms, prefix, format, options):
    "Compiler phase 5: format code"
    begin("Compiler phase 5: Formatting code")

    format.write(generated_forms, prefix, options)

    end()

def _check_options(options):
    "Check that options are valid"
    if options["optimize"]:
        warning("Optimization unavailable (will return in a future version).")
    if options["blas"]:
        warning("BLAS mode unavailable (will return in a future version).")

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
            error("Not a form: " + str(form))

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
        return QuadratureGenerator
