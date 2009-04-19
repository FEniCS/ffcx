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
__date__ = "2007-02-05 -- 2009-03-15"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.
# Modified by Dag Lindbo, 2008.
# Modified by Garth N. Wells, 2009.

__all__ = ["compile"]

# FIXME: Temporary while testing
import sys
import os

# UFL modules
from ufl.classes import Form, FiniteElementBase, Measure, Integral
from ufl.algorithms import validate_form, extract_quadrature_order, estimate_max_quadrature_order
from ufl.algorithms import extract_unique_elements, extract_basis_functions, as_form
from ufl.algorithms import *
# FFC common modules
from ffc.common.log import debug, info, warning, error, begin, end, set_level, INFO
from ffc.common.utils import product
from ffc.common.constants import UFL_OPTIONS

# FFC fem modules
from ffc.fem import create_element
from ffc.fem import create_dof_map

# FFC analysis modules
from ffc.compiler.analysis.formdata import create_ffc_coefficients
from ffc.compiler.analysis.elementdata import ElementData

# FFC form representation modules
from representation.tensor.monomialextraction import MonomialException
from representation.tensor.ufltensorrepresentation import TensorRepresentation
from representation.quadrature.uflquadraturerepresentation import QuadratureRepresentation

# FFC code generation modules
#from codegeneration.tensor import *
#from codegeneration.quadrature import *
from codegeneration.common.uflcodegenerator import generate_common_code
from codegeneration.common.integrals import generate_combined_code
from codegeneration.tensor import UFLTensorGenerator
from codegeneration.quadrature import UFLQuadratureGenerator

#from codegeneration.common.finiteelement import *
#from codegeneration.common.dofmap import *

# FFC format modules
from format.ufcformat import Format

# Form representations and code generators
if os.environ["USER"] == "logg":
    Representations = (QuadratureRepresentation, TensorRepresentation)
    CodeGenerators  = (UFLQuadratureGenerator, UFLTensorGenerator)
else:
    Representations = (QuadratureRepresentation, TensorRepresentation)
    CodeGenerators  = (UFLQuadratureGenerator, UFLTensorGenerator)

def compile(forms, prefix="Form", options=UFL_OPTIONS.copy(), global_variables=None):
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

    # Get forms and elements
    forms, elements = _extract_objects(forms)

    # Check forms
    forms = _check_forms(forms, options)

    # Check that we have at least one form
    if len(forms) == 0 and len(elements) == 0:
        info("No forms or elements specified, nothing to do.")
        return

    # Create format
    format = Format(options)

    # Compile all forms
    generated_forms = []
    for form in forms:

        # Compiler stage 1: analyze form
        form_data = analyze_form(form, options, global_variables)
        #form_datas += [form_data]

        # Compiler stage 2: compute form representations
        representations = compute_form_representations(form_data, options)

        # Compiler stage 3: optimize form representation
        optimize_form_representation(form_data)

        # Compiler stage 4: generate form code
        form_code = generate_form_code(form_data, representations, prefix, format.format, options)

        # Add to list of codes
        generated_forms += [(form_code, form_data)]

    # Generate code for elements, will only be generated if no forms were specified
    if elements:
        generated_forms += generate_element_code(elements, format.format)

    # Compiler stage 5: format code
    format_code(generated_forms, prefix, format, options)

    return
#    return (form_datas, form_representations)

def analyze_form(form, options, global_variables):
    "Compiler stage 1."
    
    begin("Compiler stage 1: Analyzing form")

    # Validate form
    validate_form(form)

    # Extract form data
    form_data = form.form_data()
    info(str(form_data))

    # Attach FFC elements and dofmaps
    form_data.ffc_elements = [create_element(element) for element in form_data.elements]
    form_data.ffc_dof_maps = [create_dof_map(element) for element in form_data.elements]

    # Attach FFC coefficients
    form_data.coefficients = create_ffc_coefficients(form_data.original_functions, global_variables)

    # FIXME: Consider adding the following to ufl.FormData
    # FIXME: Also change num_functions --> num_coefficients to match UFC
    # FIXME: Check that integrals are numbered 0, 1, 2 (not 0, 5, 6) in UFL
    form_data.num_coefficients = form_data.num_functions

    # Attach number of domains for all integral types
    form_data.num_cell_domains = max([-1] + [i.measure().domain_id() for i in form.cell_integrals()]) + 1
    form_data.num_exterior_facet_domains = max([-1] + [i.measure().domain_id() for i in form.exterior_facet_integrals()]) + 1
    form_data.num_interior_facet_domains = max([-1] + [i.measure().domain_id() for i in form.interior_facet_integrals()]) + 1

    # FIXME: Remove these, only here for compatibility with old ufcformat.py
    form_data.num_cell_integrals = form_data.num_cell_domains
    form_data.num_exterior_facet_integrals = form_data.num_exterior_facet_domains
    form_data.num_interior_facet_integrals = form_data.num_interior_facet_domains

    # Attach signature for convenience and reuse
    form_data.signature = form_data.form.signature()

    # Attach number of entries in element tensor
    dims = [create_element(v.element()).space_dimension() for v in extract_basis_functions(form)]
    dims_interior = [create_element(v.element()).space_dimension()*2 for v in extract_basis_functions(form)]
    form_data.num_entries = product(dims)
    form_data.num_entries_interior = product(dims_interior)

    # Attach number of facets
    form_data.num_facets = form_data.ffc_elements[0].num_facets()

    end()
    return form_data

def compute_form_representations(form_data, options):
    "Compiler stage 2."

    begin("Compiler stage 2: Computing form representation(s)")
    representations = [Representation(form_data) for Representation in Representations]
    end()
    
    return representations
    
def optimize_form_representation(form_data):
    "Compiler stage 3."
    
    begin("Compiler stage 3: Optimizing form representation")
    info("Optimization currently broken (to be fixed).")
    end()

def generate_form_code(form_data, representations, prefix, format, options):
    "Compiler stage 4."
    
    begin("Compiler stage 4: Generating code")

    # Generate common code like finite elements, dof map etc.
    common_code = generate_common_code(form_data, format)

    # Generate code for integrals using quadrature
    codes = []
    for (i, CodeGenerator) in enumerate(CodeGenerators):
        code_generator = CodeGenerator(options)
        codes.append(code_generator.generate_integrals(representations[i], format))

    # Loop all subdomains of integral types and combine code
    combined_code = generate_combined_code(codes, form_data, prefix, format)

    # Collect generated code
    code = {}
    code.update(common_code)
    code.update(combined_code)

    end()
    return code

def generate_element_code(elements, format):
    "Compiler stage 4."
    # Create element_data and common code
    form_data = ElementData([create_element(e) for e in elements])
    return [(generate_common_code(form_data, format), form_data)]

def format_code(generated_forms, prefix, format, options):
    "Compiler stage 5."
    
    begin("Compiler stage 5: Formatting code")
    format.write(generated_forms, prefix, options)
    end()

def _check_options(options):
    "Initial check of options."

    if "optimize" in options:
        warning("Optimization unavailable (will return in a future version).")
    if "blas" in options:
        warning("BLAS mode unavailable (will return in a future version).")
    if "quadrature_points" in options:
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
            representation = metadata["ffc_representation"]
        elif key == "quadrature_order":
            quadrature_order = metadata["quadrature_order"]
        else:
            warning("Unrecognized option '%s' for integral metadata." % key)

    # Check metadata
    valid_representations = ["tensor", "quadrature", "auto"]
    if not representation in valid_representations:
        error("Unrecognized form representation '%s', must be one of %s.",
              representation, ", ".join("'%s'" % r for r in valid_representations))
    if not ((isinstance(quadrature_order, int) and quadrature_order >= 0) or quadrature_order == "auto"):
        error("Illegal quadrature order '%s' for integral, must be a nonnegative integer or 'auto'.",
              str(quadrature_order))

    # Automatically select metadata if "auto" is selected
    if representation == "auto":
        representation = _auto_select_representation(integral)
    if quadrature_order == "auto":
        quadrature_order = _auto_select_quadrature_order(integral)

    # Create new measure with updated metadata
    metadata = {"quadrature_order": quadrature_order, "ffc_representation": representation}
    measure = integral.measure().reconstruct(metadata=metadata)

    return Integral(integral.integrand(), measure)

def _extract_objects(objects):
    "Extract forms and elements from list of objects."

    # Check that we get a list of objects
    if not isinstance(objects, (list, tuple)):
        objects = [objects]

    # Iterate over objects and extract forms and elements
    forms = []
    elements = []
    for object in objects:
        if isinstance(object, FiniteElementBase):
            elements.append(object)        
        elif not object is None:
            forms.append(as_form(object))

    # Only compile element(s) when there are no forms
    if len(forms) > 0 and len(elements) > 0:
        elements = []

    return (forms, elements)

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
        raise RuntimeError, 'Unknown form representation: "%s"' % option

def __select_code_generator(form_representation):
    "Select code generator"

    if form_representation == "tensor":
        return TensorGenerator
    else:
        return UFLQuadratureGenerator
    "Select code generator"

    if form_representation == "tensor":
        return TensorGenerator
    else:
        return UFLQuadratureGenerator

def _auto_select_representation(integral):
    "Automatically select the best representation for integral."

    # FIXME: Implement this
    info("Automatic selection of representation not implemented, defaulting to quadrature.")
    return "quadrature"

def _auto_select_quadrature_order(integral):
    "Automatically select the appropriate quadrature order for integral."

    # FIXME: to take into account Quadrature elements
    if any(e.family() == "Quadrature" for e in extract_unique_elements(integral)):
        quadrature_order = estimate_max_quadrature_order(integral._integrand)
    else:
        quadrature_order = max(estimate_max_quadrature_order(integral._integrand),\
                               estimate_quadrature_order(integral))

    return quadrature_order
