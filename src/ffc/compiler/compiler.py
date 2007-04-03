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
__date__ = "2007-02-05 -- 2007-03-23"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Kristian B. Oelgaard 2007

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
from language.reassignment import *

# FFC analysis modules
from analysis.checks import *
from analysis.formdata import *
from analysis.elementdata import *

# FFC form representation modules
from representation.tensor import *
from representation.quadrature import *

# FFC code generation modules
from codegeneration.tensor import *
from codegeneration.quadrature import *
from codegeneration.common.finiteelement import *
from codegeneration.common.dofmap import *

# FFC format modules
from format import ufcformat
from format import dolfinformat

def compile(forms, prefix = "Form", output_language = FFC_LANGUAGE, options = FFC_OPTIONS, \
            representation = FFC_REPRESENTATION):
    "Compile the given forms and/or elements for the given language"

    # Check input
    (forms, elements) = preprocess_forms(forms)
    if len(forms) == 0 and len(elements) == 0:
        debug("No forms or elements specified, nothing to do.")
        return

    # Compile forms
    if len(forms) > 0:
        compile_forms(forms, prefix, output_language, options, representation)

    # Compile elements, but only if there are no forms
    if len(elements) > 0 and len(forms) == 0:
        compile_elements(elements, prefix, output_language, options)

def compile_forms(forms, prefix = "Form", output_language = FFC_LANGUAGE, options = FFC_OPTIONS, \
                  representation = FFC_REPRESENTATION):
    "Compile the given forms for the given language"

    # Check form input
    if len(forms) == 0:
        debug("No forms specified, nothing to do.")
        return
    
    # Choose format
    format = __choose_format(output_language)

    # Iterate over forms for stages 1 - 4
    generated_forms = []
    for form in forms:

        # Compiler phase 1: analyze form
        form_data = analyze_form(form)

        # Compiler phase 2: compute form representation
        form_representation = compute_form_representation(form_data, representation)

        # Compiler phase 3: optimize form representation
        optimize_form_representation(form)

        # Compiler phase 4: generate form code
        form_code = generate_form_code(form_data, form_representation, format.format)

        # Add to list of codes
        generated_forms += [(form_code, form_data)]

    # Compiler phase 5: format code
    format_code(generated_forms, prefix, format, options)

def compile_elements(elements, prefix = "Element", output_language = FFC_LANGUAGE, options = FFC_OPTIONS):
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

    # Choose format and code generator
    format = __choose_format(output_language)
    CodeGenerator = __choose_code_generator()
    code_generator = CodeGenerator()

    # Generate code
    element_code = code_generator.generate_element_code(element_data, format.format)

    # Compiler phase 5: format code
    debug_end()
    debug_begin("Compiler phase 5: Formatting code")
    format.write([(element_code, element_data)], prefix, options)
    debug_end()
    
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

def analyze_form(form):
    "Compiler phase 1: analyze form"
    debug_begin("Phase 1: Analyzing form")

    # Check validity of form
    check_form(form)

    # Reassign form indices
    reassign_indices(form)

    # Check validity of form again
    check_form(form)

    # Extract form data
    form_data = FormData(form)

    # Print a short summary
    debug("")
    debug("Rank of form: %d" % form_data.rank)
    debug("Coefficients: %d" % form_data.num_coefficients)
    debug("Arguments:    %d" % form_data.num_arguments)
    debug("Terms:        %d" % form_data.num_terms)
    
    debug_end()
    return form_data

def compute_form_representation(form_data, representation):
    "Compiler phase 2: compute form representation"
    debug_begin("Compiler phase 2: Computing form representation")

    # Choose representation
    Representation = __choose_representation(representation)

    # Compute form representation
    form_representation = Representation(form_data)

    debug_end()
    return form_representation
    
def optimize_form_representation(form):
    "Compiler phase 3: optimize form representation"
    debug_begin("Compiler phase 3: Computing optimization")

    debug("Not implemented")

    debug_end()

def generate_form_code(form_data, form_representation, format):
    "Compiler phase 4: generate code"
    debug_begin("Compiler phase 4: Generating code")

    # Choose code generator
#    CodeGenerator = __choose_code_generator()
    if isinstance(form_representation, TensorRepresentation):
        CodeGenerator = TensorGenerator
    else:
        CodeGenerator = QuadratureGenerator

    # Generate code
    code_generator = CodeGenerator()
    code = code_generator.generate_form_code(form_data, form_representation, format)
        
    debug_end()
    return code

def format_code(generated_forms, prefix, format, options):
    "Compiler phase 5: format code"
    debug_begin("Compiler phase 5: Formatting code")

    # Format the pre-generated code
    format.write(generated_forms, prefix, options)

    debug_end()

def __choose_format(output_language):
    "Choose format from specified language"

    if output_language.lower() == "ufc":
        return ufcformat
    elif output_language.lower() == "dolfin":
        return dolfinformat
    else:
        raise RuntimeError, "Don't know how to compile code for language \"%s\"." % output_language

def __choose_representation(representation):
    "Choose form representation"

    if representation == "tensor":
        return TensorRepresentation
    else:
    # Hint: do something differently for quadrature here
        return QuadratureRepresentation

#def __choose_code_generator():
#    "Choose code generator"
    
    # Hint: do something differently for quadrature here
#    return TensorGenerator
