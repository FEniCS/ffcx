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
__date__ = "2007-02-05 -- 2009-05-04"
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
from ufl.classes import Form, FiniteElementBase, Measure, Integral, Function
from ufl.algorithms import validate_form, extract_max_quadrature_element_degree, estimate_max_polynomial_degree
from ufl.algorithms import extract_unique_elements, extract_basis_functions, as_form
from ufl.algorithms import *

# FFC common modules
from ffc.common.log import debug, info, warning, error, begin, end, set_level, INFO
from ffc.common.utils import product
from ffc.common.constants import FFC_OPTIONS

# FFC fem modules
from ffc.fem import create_element
from ffc.fem import create_dof_map

# FFC form representation modules
from representation.tensor.monomialextraction import MonomialException
from representation.tensor.tensorrepresentation import TensorRepresentation
from representation.quadrature.quadraturerepresentation import QuadratureRepresentation

# FFC code generation modules
from codegeneration.common.codegenerator import generate_common_code
from codegeneration.common.integrals import generate_combined_code
from codegeneration.tensor import TensorGenerator
from codegeneration.quadrature import QuadratureGenerator

# FFC format modules
from format.ufcformat import Format

# Form representations and code generators
Representations = (QuadratureRepresentation, TensorRepresentation)
CodeGenerators  = (QuadratureGenerator, TensorGenerator)

def compile(forms, prefix="Form", options=FFC_OPTIONS.copy(), global_variables=None):
    """This is the main interface to FFC. The input argument must be
    either a single UFL Form object or a list of UFL Form objects.
    For each form, FFC generates C++ code conforming to the UFC
    interface. The generated code is collected in a single C++ header
    file or, optionally, a pair of C++ header and implementation
    files. For detailed documentation of available options, refer to
    the FFC user manual."""

    # Check options
    options = _check_options(options)

    # Extract objects to compile
    forms, elements = _extract_objects(forms)

    # Check that we have at least one object
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

    info("Code generation complete.")
    return
#    return (form_datas, form_representations)

def analyze_form(form, options, global_variables):
    "Compiler stage 1."
    
    begin("Compiler stage 1: Analyzing form")

    # Validate form
    validate_form(form)

    # Extract form data
    form_data = form.form_data()
    form = form_data.form
    info(str(form_data))

    # Extract integral metadata
    form_data.metadata = _extract_metadata(form, options)

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
    form_data = ElementData(elements)
    return [(generate_common_code(form_data, format), form_data)]

def format_code(generated_forms, prefix, format, options):
    "Compiler stage 5."
    
    begin("Compiler stage 5: Formatting code")
    format.write(generated_forms, prefix, options)
    end()

def _check_options(options):
    "Initial check of options."

    #if "optimize" in options:
    #    warning("Optimization unavailable (will return in a future version).")
    if "blas" in options:
        warning("BLAS mode unavailable (will return in a future version).")
    if "quadrature_points" in options:
        warning("Option 'quadrature_points' has been replaced by 'quadrature_order'.")

    return options

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

def _extract_metadata(form, options):
    "Check metadata for integral and return new integral with proper metadata."

    metadata = {}

    # Iterate over integrals
    for integral in form.integrals():

        # Set default values for metadata
        representation = options["representation"]
        quadrature_order = options["quadrature_order"]

        # Get metadata for integral (if any)
        integral_metadata = integral.measure().metadata() or {}
        for (key, value) in integral_metadata.iteritems():
            if key == "ffc_representation":
                representation = integral_metadata["ffc_representation"]
            elif key == "quadrature_order":
                quadrature_order = integral_metadata["quadrature_order"]
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

        # Set metadata for integral
        metadata[integral] = {"quadrature_order": quadrature_order, "ffc_representation": representation}

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

    # Use maximum degree of quadrature element if any
    quadrature_degree = extract_max_quadrature_element_degree(integral)

    # Otherwise, estimate polynomial degree (may not be a polymomial)
    if quadrature_degree is None:
        quadrature_degree = estimate_max_polynomial_degree(integral)

    # Set quadrature order to polynomial degree
    quadrature_order = quadrature_degree

    return quadrature_order

# FIXME: Old stuff below needs to be cleaned up

def create_ffc_coefficients(ufl_functions, global_variables):
    "Try to convert UFL functions to FFC Coefficients"

    class Coefficient:
        def __init__(self, element, name):
            self.element = element
            self.name = name

    # Extract names of all coefficients
    coefficient_names = {}
    if not global_variables is None:
        for name in global_variables:
            variable = global_variables[name]
            if isinstance(variable, Function):
                coefficient_names[variable] = str(name)

    # Create Coefficients and set name
    ffc_coefficients = []
    for function in ufl_functions:
        element = function.element()
        if function in coefficient_names:
            name = coefficient_names[function]
        else:
            name = "w" + str(function.count())
        ffc_coefficients.append(Coefficient(create_element(element), name))

    return ffc_coefficients

class ElementData:
    """This class holds meta data for a list of elements. It has the
    same attributes as FormData, but only the elements and their dof
    maps are available."""

    def __init__(self, elements):
        "Create element for list of elements"

        debug("Extracting element data...")

        self.form                         = None
        self.signature                    = None
        self.rank                         = -1
        self.num_coefficients             = 0
        self.num_arguments                = len(elements)
        self.num_terms                    = 0
        self.num_cell_integrals           = 0
        self.num_exterior_facet_integrals = 0
        self.num_interior_facet_integrals = 0
        self.elements                     = [create_element(element) for element in elements]
        self.dof_maps                     = [create_dof_map(element) for element in elements]
        self.ffc_elements                 = self.elements
        self.ffc_dof_maps                 = self.dof_maps
        self.coefficients                 = []

        debug("done")
