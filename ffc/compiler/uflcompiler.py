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
__date__ = "2007-02-05 -- 2009-03-05"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.
# Modified by Dag Lindbo, 2008.

__all__ = ["compile"]

# FIXME: Temporary while testing
import sys
import os

# UFL modules
from ufl.classes import Form, FiniteElementBase, Measure, Integral
from ufl.algorithms import validate_form, extract_quadrature_order, estimate_quadrature_order
from ufl.algorithms import MonomialException, extract_unique_elements

# FFC common modules
from ffc.common.log import debug, info, warning, error, begin, end, set_level, INFO
from ffc.common.constants import FFC_OPTIONS

# FIXME: Remove this?
from analysis.formdata import FormData as FFCFormData

# FFC form representation modules
from representation.tensor.ufltensorrepresentation import TensorRepresentation
from representation.quadrature.uflquadraturerepresentation import QuadratureRepresentation

# FFC code generation modules
#from codegeneration.tensor import *
#from codegeneration.quadrature import *
from codegeneration.common.codegenerator import CodeGenerator
from codegeneration.quadrature import UFLQuadratureGenerator
#from codegeneration.common.finiteelement import *
#from codegeneration.common.dofmap import *

# FFC format modules
from format.uflufcformat import Format

def compile(forms, prefix="Form", options=FFC_OPTIONS, global_variables=None):
    """This is the main interface to FFC. The input argument must be
    either a single UFL Form object or a list of UFL Form objects.
    For each form, FFC generates C++ code conforming to the UFC
    interface. The generated code is collected in a single C++ header
    file or, optionally, a pair of C++ header and implementation
    files. For detailed documentation of available options, refer to
    the FFC user manual."""

    # Set log level
    set_level(INFO)

    # Check options
    options = _check_options(options)

    # Check forms
    forms = _check_forms(forms, options)

    # Check that we have at least one form
    if len(forms) == 0:
        info("No forms specified, nothing to do.")
        return

    # Choose format
    format = Format(options)

    # Compile all forms
    generated_forms = []
    for form in forms:

        # Compiler stage 1: analyze form
        form_data = analyze_form(form, options)
        #form_datas += [form_data]

        # Compiler stage 2: compute form representation
        tensor_representation, quadrature_representation =\
            compute_form_representation(form, options)

        # FIXME: Temporary while debugging tensor representation
        if os.environ["USER"] == "logg":
            continue

        # Compiler stage 3: optimize form representation
        # TODO: Switch this back on? I guess the argument should only be the
        # tensor_representation since it only applies to this.
        # Anders: There should be a check whether or not the form can be optimized
        # so for now it only applies to tensor representation (unless there is some
        # optimization that can be done for quadrature that is not done by default)
        #optimize_form_representation(form)

        # Compiler stage 4: generate form code
        #form_code = generate_form_code(form_data, form_representation, options["representation"], format.format)
        ffc_form_data = FFCFormData(None, global_variables, ufl_form_data=form_data)
        form_code = generate_form_code(ffc_form_data, tensor_representation,\
                          quadrature_representation, format.format)

        # Add to list of codes
        generated_forms += [(form_code, ffc_form_data)]

    if os.environ["USER"] == "logg":
        return

    # Compiler stage 5: format code
    format_code(generated_forms, prefix, format, options)

    return
#    return (form_datas, form_representations)

def analyze_form(form, options):
    "Compiler stage 1."
    
    begin("Compiler stage 1: Analyzing form")

    # Validate form
    validate_form(form)

    # Extract form data
    form_data = form.form_data()
    info(str(form_data))

    end()
    return form_data

def compute_form_representation(form, options):
    "Compiler stage 2."
    
    begin("Compiler stage 2: Computing form representation")

    # FIXME: Temporary while debugging tensor representation
    if os.environ["USER"] == "logg":

        try:
            tensor_representation = TensorRepresentation(form)
        except MonomialException, exception:
            warning("Tensor representation failed. " + exception.message)
            info("Falling back to quadrature.")
            sys.exit(1)

        return (tensor_representation, tensor_representation)

    else:
        # Compute quadrature representation
        quadrature_representation = UFLQuadratureRepresentation(form)
        tensor_representation = quadrature_representation

    end()
    return tensor_representation, quadrature_representation
    
def optimize_form_representation(form):
    "Compiler stage 3."
    
    begin("Compiler stage 3: Optimizing form representation")
    info("Optimization currently broken (to be fixed).")
    end()

def generate_form_code(form_data, tensor_representation, quadrature_representation, format):
    "Compiler stage 4."
    
    begin("Compiler stage 4: Generating code")

    # Generate common code like finite elements, dof map etc.
    common_generator = CodeGenerator()
    code = common_generator.generate_form_code(form_data, None, format, ufl_code=True)

    # We need both generators
#    tensor_generator = UFLTensorGenerator()
    tensor_generator = UFLQuadratureGenerator()
    quadrature_generator = UFLQuadratureGenerator()

    # Generate code for integrals using quadrature
    quadrature_code = quadrature_generator.generate_cell_integrals(quadrature_representation, format)
    quadrature_code.update(quadrature_generator.generate_exterior_facet_integrals(quadrature_representation, format))
    quadrature_code.update(quadrature_generator.generate_interior_facet_integrals(quadrature_representation, format))

    # Mock tensor code (remove when tensor representation generates code)
    tensor_code = {}
#    tensor_code.update(tensor_generator.generate_cell_integrals(quadrature_representation, format))
#    tensor_code.update(tensor_generator.generate_exterior_facet_integrals(quadrature_representation, format))
#    tensor_code.update(tensor_generator.generate_interior_facet_integrals(quadrature_representation, format))

    # Get any kind of code for reseting the element tensor, it just needs to be
    # generated once by the codegenerators
    reset_code = (quadrature_generator.reset_code or tensor_generator.reset_code)
    reset_code_restricted = (quadrature_generator.reset_code_restricted or tensor_generator.reset_code_restricted)

    # Loop all subdomains of integral types and combine code
    for i in range(form_data.num_cell_integrals):
        _combine_code(code, tensor_code, quadrature_code, ("cell_integral", i), reset_code, False)
    for i in range(form_data.num_exterior_facet_integrals):
        _combine_code(code, tensor_code, quadrature_code, ("exterior_facet_integral", i), reset_code, True)
    for i in range(form_data.num_interior_facet_integrals):
        _combine_code(code, tensor_code, quadrature_code, ("interior_facet_integral", i), reset_code_restricted, True)

    end()
    return code

def _combine_code(code, tensor_code, quadrature_code, key, reset_code, facet_integral):
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

def _auto_select_representation(integral):
    "Automatically select the best representation for integral."

    # FIXME: Implement this
    return "quadrature"

def _auto_select_quadrature_order(integral):
    "Automatically select the appropriate quadrature order for integral."

    # FIXME: Improve algorithms in UFL. In the mean time this is a dirty hack
    # FIXME: to take into account Quadrature elements
    if any(e.family() == "Quadrature" for e in extract_unique_elements(integral)):
        quadrature_order = extract_quadrature_order(integral)
    else:
        quadrature_order = max(extract_quadrature_order(integral),\
                               estimate_quadrature_order(integral))

    return quadrature_order

def _check_options(options):
    "Initial check of options."
    
    if options["optimize"]:
        warning("Optimization unavailable (will return in a future version).")
    if options["blas"]:
        warning("BLAS mode unavailable (will return in a future version).")
    if options["quadrature_points"]:
        warning("Option 'quadrature_points' has been replaced by 'quadrature_order'.")

    return options

def _check_forms(forms, options):
    "Initial check of forms."

    # Check that we get a list
    if not (isinstance(forms, list) or isinstance(forms, tuple)):
        forms = [forms]

    # Ignore None
    forms = [form for form in forms if not form is None]
    
    # Check that all arguments are UFL forms and ignore None
    for form in forms:
        if not isinstance(form, Form):
            error("Unable to compile, object is not a UFL form: " + str(form))

    # Check form metadata
    for (i, form) in enumerate(forms):
        new_form = Form([])
        for integral in form.integrals():
            new_form += Form([_check_metadata(integral, options)])
        forms[i] = new_form

    return forms

def _check_metadata(integral, options):
    "Check metadata for integral and return new integral with proper metadata."

    # Set default values for metadata
    representation = options["representation"]
    quadrature_order = options["quadrature_order"]

    # Get metadata for integral (if any)
    metadata = integral.measure().metadata() or {}
    for (key, value) in metadata.iteritems():
        if key == "ffc_representation":
            representation = metadata["representation"]
        elif key == "quadrature_order":
            quadrature_order = metadata["quadrature_order"]
        else:
            warning("Unrecognized option '%s' for integral metadata." % key)

    # Check metadata
    valid_representations = ["tensor", "quadrature", "automatic"]
    if not representation in valid_representations:
        error("Unrecognized form representation '%s', must be one of %s.",
              representation, ", ".join('%s' % r for r in valid_representations))
    if not ((isinstance(quadrature_order, int) and quadrature_order >= 0) or quadrature_order == "automatic"):
        error("Illegal quadrature order %s for integral, must be a nonnegative integer or 'automatic'.",
              str(quadrature_order))

    # Automatically select metadata if "automatic" is selected
    if representation == "automatic":
        representation = _auto_select_representation(integral)
    if quadrature_order == "automatic":
        quadrature_order = _auto_select_quadrature_order(integral)

    # Create new measure with updated metadata
    metadata = {"quadrature_order": quadrature_order, "ffc_representation": representation}
    measure = integral.measure().reconstruct(metadata=metadata)

    return Integral(integral.integrand(), measure)

def _extract_objects(objects):
    "Extract forms and elements from list of objects."


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
            return TensorRepresentation
        else:
            info("no\n")
            warning("Form is is not multilinear, using quadrature representation")
            return QuadratureRepresentation
        
    elif option == "quadrature":
        return QuadratureRepresentation

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
