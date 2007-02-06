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
__date__ = "2007-02-05 -- 2007-02-05"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC compiler modules
import language

from language.reassignment import *

from analysis.checks import *
from analysis.formdata import *

#import representation
#import optimization

from codegeneration.finiteelement import *
from codegeneration.dofmap import *
from codegeneration.form import *

from format import ufcformat

#import format

def compile(form, name = "Form", output_language = FFC_LANGUAGE, options = FFC_OPTIONS):
    "Compile the given form for the given language."

    # Check that we get a Form
    if not isinstance(form, language.Form):
        raise RuntimeError, "Not a form: " + str(form)

    # Phase 1: analyze form
    form_data = analyze_form(form, name)

    # Phase 2: compute form representation
    compute_representation(form)

    # Phase 3: optimize form representation
    compute_optimization(form)

    # Choose format for stages 4 and 5
    format = __choose_format(output_language)

    # Phase 4: generate code
    generate_code(form_data, format.format)

    # Phase 5: format code
    format_code(form)

def analyze_form(form, name):
    "Analyze form"

    debug("")
    debug("Compiler phase 1: Analyzing form")
    debug("--------------------------------")
    debug("")
    debug_indent()

    # Check validity of form
    check_form(form)

    # Reassign form indices
    reassign_indices(form)

    # Check validity of form again
    check_form(form)

    # Extract form data
    form_data = FormData(form, name)
    
    debug("\n" + str(form_data), 1)
    debug_indent(-1)
    
    return form_data

def compute_representation(form):
    "Compute form representation"

    debug("")
    debug("Compiler phase 2: Computing form representation")
    debug("-----------------------------------------------")
    debug_indent()

    debug("Not implemented")

    debug_indent(-1)
    
def compute_optimization(form):
    "Compute form representation"

    debug("")
    debug("Compiler phase 3: Computing optimization")
    debug("----------------------------------------")
    debug_indent()

    debug("Not implemented")

    debug_indent(-1)

def generate_code(form_data, format):
    "Generate code"

    debug("")
    debug("Compiler phase 4: Generating code")
    debug("---------------------------------")
    debug_indent()

    code = {}

    # Set name
    code["name"] = form_data.name

    # Set number of arguments
    #code["num_arguments"] = form_data.num_arguments

    # Generate code for finite elements
    for i in range(len(form_data.elements)):
        code[("finite_element", i)] = generate_finite_element(form_data.elements[i], format)

    # Generate code for dof maps
    for i in range(len(form_data.dof_maps)):
        code[("dof_map", i)] = generate_dof_map(form_data.dof_maps[i], format)

    # Generate code for form
    code["form"] = generate_form(form_data, format)

    print code
    
    debug_indent(-1)

def format_code(form):
    "Format code"

    debug("")
    debug("Compiler phase 5: Formatting code")
    debug("---------------------------------")
    debug_indent()

    debug("Not implemented")

    debug_indent(-1)

def __choose_format(output_language):
    "Choose format from specified language."

    if output_language.lower() == "ufc":
        return ufcformat
    else:
        raise RuntimeError, "Don't know how to compile code for language \"%s\".", output_language
