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
__date__ = "2007-02-05 -- 2007-02-07"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC compiler modules
import language

from language.reassignment import *

# FFC analysis modules
from analysis.checks import *
from analysis.formdata import *

from representation.tensor import *

#import representation
#import optimization

# FFC code generation modules
from codegeneration.finiteelement import *
from codegeneration.dofmap import *
from codegeneration.cellintegral import *
from codegeneration.exteriorfacetintegral import *
from codegeneration.interiorfacetintegral import *
from codegeneration.form import *

# FFC format modules
from format import ufcformat

def compile(form, name = "Form", output_language = FFC_LANGUAGE, options = FFC_OPTIONS):
    "Compile the given form for the given language."

    # Check that we get a Form
    if not isinstance(form, language.Form):
        raise RuntimeError, "Not a form: " + str(form)

    # Compiler phase 1: analyze form
    form_data = analyze_form(form, name)

    # Compiler phase 2: compute form representation
    form_representation = compute_representation(form_data)

    # Compiler phase 3: optimize form representation
    compute_optimization(form)

    # Choose format for stages 4 and 5
    format = __choose_format(output_language)

    # Compiler phase 4: generate code
    code = generate_code(form_data, form_representation, format.format)

    # Compiler phase 5: format code
    format_code(code, form_data, format, options)

def analyze_form(form, name):
    "Compiler phase 1: analyze form"
    debug_begin("Phase 1: Analyzing form")

    # Check validity of form
    check_form(form)

    # Reassign form indices
    reassign_indices(form)

    # Check validity of form again
    check_form(form)

    # Extract form data
    form_data = FormData(form, name)

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

    # At this point, we need to choose the type of representation, but
    # currently only the tensor representation is implemented.
    # Hint: do something differently for quadrature here.

    # Compute tensor representation
    form_representation = TensorRepresentation(form_data)

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

    code = {}

    # Generate code for finite elements
    debug("Generating code for finite elements...")
    for i in range(len(form_data.elements)):
        code[("finite_element", i)] = generate_finite_element(form_data.elements[i], format)
    debug("done")

    # Generate code for dof maps
    debug("Generating code for dof maps...")
    for i in range(len(form_data.dof_maps)):
        code[("dof_map", i)] = generate_dof_map(form_data.dof_maps[i], format)
    debug("done")

    # FIXME: Should be multiple integrals for each type

    # Generate code for cell integral
    debug("Generating code for cell integral...")
    code["cell_integral"] = generate_cell_integral(form_representation, format)
    debug("done")

    # Generate code for cell exterior facet integral
    debug("Generating code for exterior facet integral...")
    code["exterior_facet_integral"] = generate_exterior_facet_integral(form_representation, format)
    debug("done")

    # Generate code for cell interior facet integral
    debug("Generating code for interior facet integral...")
    code["interior_facet_integral"] = generate_interior_facet_integral(form_representation, format)
    debug("done")

    # Generate code for form
    debug("Generating code for form...")
    code["form"] = generate_form(form_data, format)
    debug("done")

    debug_end()
    return code

def format_code(code, form_data, format, options):
    "Compiler phase 5: format code"
    debug_begin("Compiler phase 5: Formatting code")

    # Format the pre-generated code
    format.write(code, form_data, options)

    debug_end()

def __choose_format(output_language):
    "Choose format from specified language."

    if output_language.lower() == "ufc":
        return ufcformat
    else:
        raise RuntimeError, "Don't know how to compile code for language \"%s\".", output_language
