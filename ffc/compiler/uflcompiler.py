"""This is the compiler, acting as the main interface for compilation
of forms and breaking the compilation into several sequential phases,
each represented by a separate module:

   0. language        -  expressing the form in the form language
   1. analysis        -  simplifying and preprocessing the form
   2. representation  -  computing a representation of the form
   3. optimization    -  optimizing the form representation
   4. codegeneration  -  generating code according to a format
   5. format          -  writing the generated code to file

"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-05 -- 2008-10-21"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2007
# Modified by Dag Lindbo, 2008

# UFL modules
from ufl.classes import Form, FiniteElementBase
from ufl.algorithms import FormData, is_multilinear

# FFC common modules
from ffc.common.debug import debug, warning, debug_begin, debug_end
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
    "Compile the given forms and/or elements"

    warning("UFL compiler is experimental. In particular, it doesn't work yet.")

    # Check options
    _check_options(options)

    # Extract forms and elements
    (forms, elements) = _extract_objects(objects)
    if len(forms) == 0 and len(elements) == 0:
        debug("No forms or elements specified, nothing to do.")
        return

    # Compile forms
    (form_data, form_representation) = _compile_forms(forms, prefix, options)

    # Compile elements
    #_compile_elements(elements, prefix, options)

    return (form_data, form_representation)

def _compile_forms(forms, prefix, options):
    "Compile the given forms"

    # Check input
    if len(forms) == 0:
        return (None, None)

    # Choose format
    format = _choose_format(options["language"])
    format.init(options)

    # Iterate over forms for stages 1 - 4
    generated_forms = []
    form_datas = []
    form_representations = []
    for form in forms:

        # Compiler phase 1: analyze form
        print "1"
        form_data = analyze_form(form)
        form_datas += [form_data]

        # Compiler phase 2: compute form representation
        form_representation = compute_form_representation(form_data, options)
        form_representations += [form_representation]

        # Compiler phase 3: optimize form representation
        #optimize_form_representation(form)

        # Compiler phase 4: generate form code
        #form_code = generate_form_code(form_data, form_representation, options["representation"], format.format)

        # Add to list of codes
        #generated_forms += [(form_code, form_data)]

    # Compiler phase 5: format code
    #_format_code(generated_forms, prefix, format, options)

    return (form_datas, form_representations)

def __compile_elements(elements, prefix="Element", options=FFC_OPTIONS):
    "Compile the given elements for the given language"

    # Check input
    if len(forms) == 0:
        return

    # Compiler phase 1: analyze form
    debug_begin("Compiler phase 1: Analyzing elements")
    element_data = ElementData(elements)
    debug_end()

    # Go directly to phase 4, code generation
    debug_begin("Compiler phase 2-3: Nothing to do")
    debug("-")
    debug_end()
    debug_begin("Compiler phase 4: Generating code")

    # Choose format
    format = __choose_format(options["language"])
    format.init(options)

    # Choose code generator
    CodeGenerator = __choose_code_generator(options["representation"])
    code_generator = CodeGenerator()

    # Generate code
    element_code = code_generator.generate_element_code(element_data, format.format)

    # Compiler phase 5: format code
    debug_end()
    debug_begin("Compiler phase 5: Formatting code")
    format.write([(element_code, element_data)], prefix, options)
    debug_end()
    
def analyze_form(form):
    "Compiler phase 1: analyze form"
    debug_begin("Phase 1: Analyzing form")

    # Analyze form and extract form data
    form_data = FormData(form)
    debug(str(form_data))

    debug_end()
    return form_data

def compute_form_representation(form_data, options):
    "Compiler phase 2: compute form representation"
    debug_begin("Compiler phase 2: Computing form representation")

    # Choose representation
    Representation = _choose_representation(form_data.form, options)

    # Compute form representation
    form_representation = Representation(form_data, int(options["quadrature_points"]))

    debug_end()
    return form_representation
    
def optimize_form_representation(form):
    "Compiler phase 3: optimize form representation"
    debug_begin("Compiler phase 3: Computing optimization")

    debug("Not implemented")

    debug_end()

def generate_form_code(form_data, form_representation, representation, format):
    "Compiler phase 4: generate code"
    debug_begin("Compiler phase 4: Generating code")

    # Choose code generator
    CodeGenerator = __choose_code_generator(representation)

    # Generate code
    code_generator = CodeGenerator()
    code = code_generator.generate_form_code(form_data, form_representation, format)
        
    debug_end()
    return code

def format_code(generated_forms, prefix, format, options):
    "Compiler phase 5: format code"
    debug_begin("Compiler phase 5: Formatting code")

    format.write(generated_forms, prefix, options)

    debug_end()

def _check_options(options):
    "Check that options are valid"
    # FIXME: We could do more tests here
    if options["optimize"]:
        debug("*** Warning: " + "Optimization unavailable (will return in a future version)")
    if options["blas"]:
        debug("*** Warning: " + "BLAS mode unavailable (will return in a future version)")

def _extract_objects(objects):
    "Extract forms and elements from list of objects"

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
            raise RuntimeError, "Not a form: " + str(form)

    # Only compile element(s) when there are no forms
    if len(forms) > 0:
        elements = []

    return (forms, elements)

def _choose_format(language):
    "Choose format from specified language"

    # FIXME: Make format a class (since we call init())

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
        debug("Checking if form is multilinear...")
        if is_multilinear(form):
            debug("yes\n")
            return UFLTensorRepresentation
        else:
            debug("no\n")
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
