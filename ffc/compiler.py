"""
This is the compiler, acting as the main interface for compilation
of forms and breaking the compilation into several sequential stages.
The output of each stage is the input of the next stage.

Stage 0: Language, parsing
--------------------------

  Input:  Python code or .ufl file
  Output: UFL form

  This stage consists of parsing and expressing a form in the
  UFL form language.

  This stage is completely handled by UFL.

Stage 1: Analysis
-----------------

  Input:  UFL form
  Output: Preprocessed UFL form and FormData (metadata)

  This stage preprocesses the UFL form and extracts form metadata.
  It may also perform simplifications on the form.

Stage 2: Code representation
----------------------------

  Input:  Preprocessed UFL form and FormData (metadata)
  Output: Intermediate Representation (IR)

  This stage examines the input and generates all data needed for code
  generation. This includes generation of finite element basis
  functions, extraction of data for mapping of degrees of freedom and
  possible precomputation of integrals.

  Most of the complexity of compilation is handled in this stage.

  The IR is stored as a dictionary, mapping names of UFC functions to
  data needed for generation of the corresponding code.

Stage 3: Code optimization
--------------------------

  Input:  Intermediate Representation (IR)
  Output: Optimized Intermediate Representation (OIR)

  This stage examines the IR and performs optimizations.

Stage 4: Code generation
------------------------

  Input:  Optimized Intermediate Representation (OIR)
  Output: C++ code

  This stage examines the OIR and generates the actual C++ code for
  the body of each UFC function.

  The code is stored as a dictionary, mapping names of UFC functions
  to strings containing the C++ code of the body of each function.

Stage 5: Code formatting
------------------------

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
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.
# Modified by Dag Lindbo, 2008.
# Modified by Garth N. Wells, 2009.
# Last changed: 2009-12-09

__all__ = ["compile_form", "compile_element"]

# UFL modules.
from ufl.common import istr
from ufl.algorithms import preprocess
from ufl.algorithms import FormData
from ufl.algorithms import estimate_max_polynomial_degree
from ufl.algorithms import estimate_total_polynomial_degree

# FFC modules.
from log import log, begin, end, debug, info, warning, error, ffc_assert
from constants import FFC_OPTIONS
from quadratureelement import default_quadrature_degree

# FFC representation modules
from representation import compute_form_ir
from representation import compute_element_ir
from representation import compute_dofmap_ir

# FFC code generation modules
from codegeneration import generate_element_code
from codegeneration import generate_dofmap_code

# FFC formatting modules
from formatting import format_ufc

# Representation methods
methods = ("quadrature", "tensor")

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
    format_code(codes, options)

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

def analyze_form(form, options):
    "Compiler stage 1."

    begin("Compiler stage 1: Analyzing form")

    # Preprocess form
    if not form.is_preprocessed():
        form = preprocess(form)

    # Compute form data
    form_data = FormData(form)
    info(str(form_data))

    # Adjust cell and degree for elements when unspecified
    _adjust_elements(form_data)

    # Extract integral metadata
    form_data.metadata = _extract_metadata(form, options, form_data.elements)

    end()
    return form, form_data

def compute_ir(form, form_data, options):
    """Compiler stage 2.

    This stage generates an intermediate representation of the given
    form and all its elements and corresponding dofmaps. Note that
    a list of representations is generated by iterating over possible
    form representations (quadrature and tensor). For elements and
    dofmaps we always obtain a list of representations since each
    form may be associated with multiple elements and dofmaps."""

    begin("Compiler stage 2: Computing intermediate representation")

    # Generate representations for form
    begin("Computing form representations")
    ir_forms = [compute_form_ir(form, form_data, m) for m in methods]
    end()

    # Generate representations for finite elements
    begin("Computing element representations")
    ir_elements = [compute_element_ir(e) for e in form_data.elements]
    end()

    # Generate representations for dofmaps
    begin("Computing dofmap representations")
    ir_dofmaps = [compute_dofmap_ir(e) for e in form_data.elements]
    end()

    end()

    return ir_forms, ir_elements, ir_dofmaps

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
    ir_forms, ir_elements, ir_dofmaps = ir

    # Generate code for forms
    # FIXME: Not implemented

    # Generate code for elements
    code_elements = [generate_element_code(ir) for ir in ir_elements]

    # Generate code for dofmaps
    code_dofmaps = [generate_dofmap_code(ir) for ir in ir_dofmaps]

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

#def generate_element_code(elements, format):
#    "Compiler stage 4."
#    # Create element_data and common code
#    form_data = ElementData(elements)
#    return [(generate_common_code(form_data, format), form_data)]

def format_code(codes, options):
    "Compiler stage 5."

    # FIXME: Needs to be updated

    begin("Compiler stage 5: Formatting code")
    #format.write(codes, options)
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
        warning("Option 'quadrature_points' has been replaced by 'quadrature_order'.")

    return options

def _extract_metadata(form, options, elements):
    "Check metadata for integral and return new integral with proper metadata."

    metadata = {}

    # Iterate over integrals
    for integral in form.integrals():

        # Set default values for metadata
        representation = options["representation"]
        quadrature_order = options["quadrature_order"]
        quadrature_rule = options["quadrature_rule"]

        if quadrature_rule is None:
            info("Quadrature rule: default")
        else:
            info("Quadrature rule: " + str(quadrature_rule))
        info("Quadrature order: " + str(quadrature_order))

        # Get metadata for integral (if any)
        integral_metadata = integral.measure().metadata() or {}
        for (key, value) in integral_metadata.iteritems():
            if key == "ffc_representation":
                representation = integral_metadata["ffc_representation"]
            elif key == "quadrature_order":
                quadrature_order = integral_metadata["quadrature_order"]
            elif key == "quadrature_rule":
                quadrature_rule = integral_metadata["quadrature_rule"]
            else:
                warning("Unrecognized option '%s' for integral metadata." % key)

        # Check metadata
        valid_representations = ["tensor", "quadrature", "auto"]
        if not representation in valid_representations:
            error("Unrecognized form representation '%s', must be one of %s.",
                  representation, ", ".join("'%s'" % r for r in valid_representations))
        if quadrature_order != "auto":
            try:
                quadrature_order = int(quadrature_order)
                if not quadrature_order >= 0:
                    error("Illegal quadrature order '%s' for integral, must be a nonnegative integer.",
                        str(quadrature_order))
            except:
                error("Illegal quadrature order '%s' for integral, must be a nonnegative integer or 'auto'.",
                    str(quadrature_order))

        # FIXME: Change from quadrature_order --> quadrature_degree

        # Automatically select metadata if "auto" is selected
        if representation == "auto":
            representation = _auto_select_representation(integral)
        if quadrature_order == "auto":
            quadrature_order = _auto_select_quadrature_degree(integral, representation, elements)
        log(30, "Integral quadrature degree is %d." % quadrature_order)

        # No quadrature rules have been implemented yet
        if quadrature_rule:
            warning("No quadrature rules have been implemented yet, using the default from FIAT.")

        # Set metadata for integral
        metadata[integral] = {"quadrature_order": quadrature_order,
                              "ffc_representation": representation,
                              "quadrature_rule":quadrature_rule}

    return metadata

def _select_representation(form, options):
    "Select form representation"

    option = options["representation"]
    if option == "tensor":

        # Check if form is multilinear
        info("Checking if form is multilinear...")
        if is_multilinear(form):
            info("yes\n")
            return TensorRepresentation
        else:
            info("no\n")
            warning("Form is is not multilinear, using quadrature representation")
            return QuadratureRepresentation

    elif option == "quadrature":
        return QuadratureRepresentation

    else:
        error('Unknown form representation: "%s"' % option)

def _auto_select_representation(integral):
    "Automatically select the best representation for integral."

    # FIXME: Implement this
    info("Automatic selection of representation not implemented, defaulting to quadrature.")
    return "quadrature"

def _auto_select_quadrature_degree(integral, representation, elements):
    "Automatically select the appropriate quadrature degree for integral."

    # Estimate total degree of integrand
    degree = estimate_total_polynomial_degree(integral, default_quadrature_degree)

    # Use maximum quadrature element degree if any for quadrature representation
    if representation == "quadrature":
        #quadrature_elements = [e for e in elements if e.family() == "Quadrature"]
        #degree = max([degree] + [e.degree() for e in quadrature_elements])
        quadrature_degrees = [e.degree() for e in elements if e.family() == "Quadrature"]
        if quadrature_degrees != []:
            ffc_assert(min(quadrature_degrees) == max(quadrature_degrees), \
                       "All QuadratureElements in an integrand must have the same degree: %s" \
                       % str(quadrature_degrees))
            degree = quadrature_degrees[0]

    return degree

def _adjust_elements(form_data):
    "Adjust cell and degree for elements when unspecified"

    # Extract common cell
    common_cell = form_data.cell
    if common_cell.domain() is None:
        error("Missing cell definition in form.")

    # Extract common degree
    common_degree = max([element.degree() for element in form_data.elements])
    if common_degree is None:
        common_degree = default_quadrature_degree

    # Set cell and degree if missing
    for element in form_data.elements:

        # Check if cell and degree need to be adjusted
        cell = element.cell()
        degree = element.degree()
        if degree is None:
            #info("Adjusting element degree from %s to %d" % (istr(degree), common_degree))
            log(30, "Adjusting element degree from %s to %d" % (istr(degree), common_degree))
            element.set_degree(common_degree)
        if cell.domain() is None:
            #info("Adjusting element cell from %s to %s." % (istr(cell), str(common_cell)))
            log(30, "Adjusting element cell from %s to %s." % (istr(cell), str(common_cell)))
            element.set_cell(common_cell)
