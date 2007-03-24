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

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC FEM modules
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

# FFC form representation modules
from representation.tensor import *
from representation.quadrature import *

# FFC code generation modules
from codegeneration.tensor import *
from codegeneration.quadrature import *

# FFC format modules
from format import ufcformat
from format import dolfinformat

def compile(forms, prefix = "Form", output_language = FFC_LANGUAGE, options = FFC_OPTIONS):
    "Compile the given form for the given language."

    # Check form input
    forms = preprocess_forms(forms)
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
        form_representation = compute_representation(form_data)

        # Compiler phase 3: optimize form representation
        compute_optimization(form)

        # Compiler phase 4: generate code
        form_code = generate_code(form_data, form_representation, format.format)

        # Add to list of codes
        generated_forms += [(form_code, form_data)]

    # Compiler phase 5: format code
    format_code(generated_forms, prefix, format, options)

def preprocess_forms(forms):
    "Check and possibly convert form input to a list of Forms"

    # Check that we get a list of forms
    if not isinstance(forms, list):
        forms = [forms]

    # Check each form
    preprocessed_forms = []
    for form in forms:
        if isinstance(form, Form):
            preprocessed_forms += [form]
        elif isinstance(form, Monomial):
            preprocessed_forms += [Form(form)]
        elif not form == None:
            raise RuntimeError, "Not a form: " + str(form)

    return preprocessed_forms

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

def compute_representation(form_data):
    "Compiler phase 2: compute form representation"
    debug_begin("Compiler phase 2: Computing form representation")

    # Choose representation
    Representation = __choose_representation()

    # Compute form representation
    form_representation = Representation(form_data)

    debug_end()
    return form_representation
    
def compute_optimization(form):
    "Compiler phase 3: optimize form representation"
    debug_begin("Compiler phase 3: Computing optimization")

    debug("Not implemented")

    debug_end()

def generate_code(form_data, form_representation, format):
    "Compiler phase 4: generate code"
    debug_begin("Compiler phase 4: Generating code")

    # Choose code generator
    CodeGenerator = __choose_code_generator()

    # Generate code
    code_generator = CodeGenerator(form_data, form_representation, format)
    code = code_generator.generate_code()
        
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

def __choose_representation():
    "Choose form representation"
    
    # Hint: do something differently for quadrature here
    return TensorRepresentation

def __choose_code_generator():
    "Choose code generator"
    
    # Hint: do something differently for quadrature here
    return TensorGenerator
