"""
This is the compiler, acting as the main interface for compilation
of forms and breaking the compilation into several sequential stages.
The output of each stage is the input of the next stage.

Compiler stage 0: Language, parsing
-----------------------------------

  Input:  Python code or .ufl file
  Output: UFL form

  This stage consists of parsing and expressing a form in the
  UFL form language.

  This stage is completely handled by UFL.

Compiler stage 1: Analysis
--------------------------

  Input:  UFL form
  Output: Preprocessed UFL form and FormData (metadata)

  This stage preprocesses the UFL form and extracts form metadata.
  It may also perform simplifications on the form.

Compiler stage 2: Code representation
-------------------------------------

  Input:  Preprocessed UFL form and FormData (metadata)
  Output: Intermediate Representation (IR)

  This stage examines the input and generates all data needed for code
  generation. This includes generation of finite element basis
  functions, extraction of data for mapping of degrees of freedom and
  possible precomputation of integrals.

  Most of the complexity of compilation is handled in this stage.

  The IR is stored as a dictionary, mapping names of UFC functions to
  data needed for generation of the corresponding code.

Compiler stage 3: Code optimization
-----------------------------------

  Input:  Intermediate Representation (IR)
  Output: Optimized Intermediate Representation (OIR)

  This stage examines the IR and performs optimizations.

  Optimization is currently disabled as a separate stage
  but is implemented as part of the code generation for
  quadrature representation.

Compiler stage 4: Code generation
---------------------------------

  Input:  Optimized Intermediate Representation (OIR)
  Output: C++ code

  This stage examines the OIR and generates the actual C++ code for
  the body of each UFC function.

  The code is stored as a dictionary, mapping names of UFC functions
  to strings containing the C++ code of the body of each function.

Compiler stage 5: Code formatting
---------------------------------

  Input:  C++ code
  Output: C++ code files

  This stage examines the generated C++ code and formats it according
  to the UFC format, generating as output one or more .h/.cpp files
  conforming to the UFC format.

The main interface is defined by the following two functions:

  compile_form
  compile_element

The compiler stages are implemented by the following functions:

  analyze_form  (stage 1)
  compute_ir    (stage 2)
  optimize_ir   (stage 3)
  generate_code (stage 4)
  format_code   (stage 5)
"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-05"
__copyright__ = "Copyright (C) 2007-2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.
# Modified by Dag Lindbo, 2008.
# Modified by Garth N. Wells, 2009.
# Last changed: 2009-12-09

__all__ = ["compile_form", "compile_element"]

# FFC modules
from ffc.log import log, begin, end, debug, info, warning, error, ffc_assert
from ffc.constants import FFC_OPTIONS

# FFC analysis modules
from ffc.analysis import analyze_form

# FFC representation modules
from ffc.representation import compute_form_ir
from ffc.representation import compute_element_ir
from ffc.representation import compute_dofmap_ir

# FFC code generation modules
from ffc.codegeneration import generate_form_code
from ffc.codegeneration import generate_element_code
from ffc.codegeneration import generate_dofmap_code

# FFC formatting modules
from ffc.formatting import format_ufc

# Representation types
representation_types = ("quadrature", "tensor")

def compile_form(forms, prefix="Form", options=FFC_OPTIONS.copy()):
    """This function generates UFC code for a given UFL form or list
    of UFL forms."""

    # Check input arguments
    forms = _check_forms(forms)
    options = _check_options(options)
    if not forms: return

    # Storage for generated data
    form_and_data = []
    codes = []

    # Enter compiler stages 1-4 for each form
    for form in forms:

        # Stage 1: analysis
        preprocessed_form, form_data = analyze_form(form, options)

        # Stage 2: intermediate representation
        ir = compute_ir(preprocessed_form, form_data, options)

        # Stage 3: optimization
        oir = optimize_ir(ir, options)

        # Stage 4: code generation
        code = generate_code(oir, options)

        # Store data
        form_and_data.append((preprocessed_form, form_data))
        codes.append(code)

    # Stage 5: format code
    format_code(codes, prefix, options)

    info("Code generation complete.")
    return form_and_data

def compile_element(elements, prefix="Element", options=FFC_OPTIONS.copy()):
    """This function generates UFC code for a given UFL element or
    list of UFL elements."""

    # Check options
    options = _check_options(options)

    # Check that we get a list of elements
    if not isinstance(elements, (list, tuple)):
        elements = [elements]

    # Check that we have at least one element
    if len(elements) == 0:
        info("No elements specified, nothing to do.")
        return

    # Create format
    format = Format(options)

    # Compiler stage 4: generate element code
    #generated_elements = generate_element_code(elements, format.format)

    # Compiler stage 5: format code
    format_code(generated_elements, prefix, format, options)

    info("Code generation complete.")

def compute_ir(form, form_data, options):
    """Compiler stage 2.

    This stage generates an intermediate representation of the given
    form and all its elements and corresponding dofmaps. Note that
    a list of representations is generated by iterating over possible
    form representations (quadrature and tensor). For elements and
    dofmaps we always obtain a list of representations since each
    form may be associated with multiple elements and dofmaps."""

    begin("Compiler stage 2: Computing intermediate representation")

    # Compute form representation
    begin("Computing form representation")
    ir_form = compute_form_ir(form, form_data)
    end()

    # Compute representations for finite elements (unique sub elements)
    begin("Computing element representations")
    ir_elements = [compute_element_ir(e) for e in form_data.unique_sub_elements]
    end()

    # Compute representations for dofmaps (unique sub dofmaps)
    begin("Computing dofmap representations")
    ir_dofmaps = [compute_dofmap_ir(e) for e in form_data.unique_sub_elements]
    end()

    end()

    return ir_form, ir_elements, ir_dofmaps

def optimize_ir(ir, form_data):
    "Compiler stage 3."

    begin("Compiler stage 3: Optimizing intermediate representation")
    info("Optimization is currently missing")
    end()

    return ir

def generate_code(ir, options):
    "Compiler stage 4."

    begin("Compiler stage 4: Generating code")

    # Extract intermediate representations
    ir_form, ir_elements, ir_dofmaps = ir

    # Generate code for forms
    code_form = generate_form_code(ir_form, options)

    # Generate code for elements
    code_elements = [generate_element_code(ir, options) for ir in ir_elements]

    # Generate code for dofmaps
    code_dofmaps = [generate_dofmap_code(ir, options) for ir in ir_dofmaps]

    return code_form, code_elements, code_dofmaps

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

def format_code(codes, prefix, options):
    "Compiler stage 5."

    # FIXME: Needs to be updated

    begin("Compiler stage 5: Formatting code")
    format_ufc(codes, prefix, options)
    end()

def _check_forms(forms):
    "Initial check of forms."

    if not isinstance(forms, (list, tuple)):
        forms = (forms,)

    return forms

def _check_options(options):
    "Initial check of options."

    if "blas" in options:
        warning("BLAS mode unavailable (will return in a future version).")
    if "quadrature_points" in options:
        warning("Option 'quadrature_points' has been replaced by 'quadrature_degree'.")

    return options
