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
__date__ = "2007-02-05 -- 2008-09-04"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard 2007
# Modified by Dag Lindbo, 2008

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC fem modules
from ffc.fem.finiteelement import *
from ffc.fem.vectorelement import *
from ffc.fem.mixedelement import *
from ffc.fem.mixedfunctions import *

# FFC language modules
from language.algebra import *

# FFC analysis modules
from analysis.analyze import *

# FFC form representation modules
from representation.tensor.tensorrepresentation import *
from representation.quadrature.quadraturerepresentation import *

# FFC code generation modules
from codegeneration.tensor import *
from codegeneration.quadrature import *
from codegeneration.common.finiteelement import *
from codegeneration.common.dofmap import *

# FFC format modules
from format.ufcformat import Format

def compile(forms, prefix="Form", options=FFC_OPTIONS, global_variables=None):
    "Compile the given forms and/or elements"

    # Check options
    check_options(options)
    
    # Check input
    (forms, elements) = preprocess_forms(forms)
    if len(forms) == 0 and len(elements) == 0:
        debug("No forms or elements specified, nothing to do.")
        debug_end()
        return

    # Compile forms
    form_data = None
    form_representation = None
    if len(forms) > 0:
        (form_data, form_representation) = __compile_forms(forms, prefix, options, global_variables)

    # Compile elements, but only if there are no forms
    if len(elements) > 0 and len(forms) == 0:
        r1 = __compile_elements(elements, prefix, options)

    return (form_data, form_representation)

def __compile_forms(forms, prefix, options, global_variables):
    "Compile the given forms"

    # Check form input
    if len(forms) == 0:
        debug("No forms specified, nothing to do.")
        return

    # Initialise format
    format = Format(options)

    # Iterate over forms for stages 1 - 4
    generated_forms = []
    form_datas = []
    form_representations = []
    for form in forms:

        # Compiler phase 1: analyze form
        form_data = analyze_form(form, global_variables)
        form_datas += [form_data]

        # Compiler phase 2: compute form representation
        form_representation = compute_form_representation(form_data, options)
        form_representations += [form_representation]

        # Compiler phase 3: optimize form representation
        optimize_form_representation(form)

        # Compiler phase 4: generate form code
        form_code = generate_form_code(form_data, form_representation, options["representation"], format.format)

        # Add to list of codes
        generated_forms += [(form_code, form_data)]

    # Compiler phase 5: format code
    format_code(generated_forms, prefix, format, options)

    return (form_datas, form_representations)

def __compile_elements(elements, prefix, options):
    "Compile the given elements for the given language"

    # Check element input
    if len(elements) == 0:
        debug("No elements specified, nothing to do.")
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

    # Initialise format
    format = Format(options)

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
    
def check_options(options):
    "Check that options are valid"
    # FIXME: We could do more tests here
    if options["optimize"]:
        debug("*** Warning: " + "Optimization unavailable (will return in a future version)")
    if options["blas"]:
        debug("*** Warning: " + "BLAS mode unavailable (will return in a future version)")
    if options["quadrature_order"] != "automatic":
        debug("*** Warning: " + "The option 'quadrature_order' is only available for UFL forms, will be ignored")

def preprocess_forms(forms):
    "Check and possibly convert form input to a list of Forms and a list if elements"

    # Check that we get a list
    if not isinstance(forms, list):
        forms = [forms]

    # Check each form
    preprocessed_forms = []
    preprocessed_elements = []
    for form in forms:
        if isinstance(form, Form):
            preprocessed_forms += [form]
        elif isinstance(form, Monomial):
            preprocessed_forms += [Form(form)]
        elif isinstance(form, FiniteElement) or isinstance(form, MixedElement):
            preprocessed_elements += [form]
        elif not form == None:
            raise RuntimeError, "Not a form: " + str(form)

    return (preprocessed_forms, preprocessed_elements)

def analyze_form(form, global_variables):
    "Compiler phase 1: analyze form"
    debug_begin("Compiler phase 1: Analyzing form")

    # Analyze form and extract form data
    form_data = analyze(form, global_variables=global_variables)

    debug_end()
    return form_data

def compute_form_representation(form_data, options):
    "Compiler phase 2: compute form representation"
    debug_begin("Compiler phase 2: Computing form representation")

    # Choose representation
    Representation = __choose_representation(options["representation"])

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

def __choose_representation(form_representation):
    "Choose form representation"

    if form_representation == "tensor":
        return TensorRepresentation
    else:
        return QuadratureRepresentation

def __choose_code_generator(form_representation):
    "Choose code generator"

    if form_representation == "tensor":
        return TensorGenerator
    else:
        return QuadratureGenerator
