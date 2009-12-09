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
__date__ = "2007-02-05"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009.
# Modified by Dag Lindbo, 2008.
# Modified by Garth N. Wells, 2009.
# Last changed: 2009-12-08

__all__ = ["compile"]

# UFL modules
from ufl.algorithms import preprocess, FormData
from ufl.algorithms import estimate_max_polynomial_degree, estimate_total_polynomial_degree
from ufl.classes import FiniteElementBase, Integral
from ufl.common import istr

# FFC common modules
from ffc.common.log import begin, end, debug, info, warning, error, log, ffc_assert
from ffc.common.utils import product
from ffc.common.constants import FFC_OPTIONS

# FFC fem modules
from ffc.fem import create_element, create_dof_map
from ffc.fem.quadratureelement import default_quadrature_degree

# FFC form representation modules
from tensor.tensorrepresentation import TensorRepresentation
from quadrature.quadraturerepresentation import QuadratureRepresentation

## FFC code generation modules
from codegenerator import generate_common_code
from integrals import generate_combined_code
from tensor.tensorgenerator import TensorGenerator
from quadrature.quadraturegenerator import QuadratureGenerator


# Form representations and code generators
Representations = (QuadratureRepresentation, TensorRepresentation)
CodeGenerators  = (QuadratureGenerator, TensorGenerator)

def compile(objects, prefix="Form", options=FFC_OPTIONS.copy()):
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
    forms, elements = _extract_objects(objects)

    # Check that we have at least one object
    if len(forms) == 0 and len(elements) == 0:
        info("No forms or elements specified, nothing to do.")
        return

    # Create format
    format = Format(options)

    # Compile all forms
    form_and_data = []
    generated_forms = []
    for form in forms:

        # Compiler stage 1: analyze form
        form, form_data = analyze_form(form, options)
        form_and_data.append((form, form_data))

        # Compiler stage 2: compute form representations
        representations = compute_form_representations(form, form_data, options)

        # Compiler stage 3: optimize form representation
        optimize_form_representation(form, form_data)

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
    return form_and_data

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

    # Attach FFC elements and dofmaps
    form_data.ffc_elements = [create_element(element) for element in form_data.elements]
    form_data.ffc_dof_maps = [create_dof_map(element) for element in form_data.elements]

    # Attach FFC coefficients
    form_data.coefficients = create_ffc_coefficients(form_data.coefficients, form_data.coefficient_names)

    # Attach number of domains for all integral types
    form_data.num_cell_domains = max([-1] + [i.measure().domain_id() for i in form.cell_integrals()]) + 1
    form_data.num_exterior_facet_domains = max([-1] + [i.measure().domain_id() for i in form.exterior_facet_integrals()]) + 1
    form_data.num_interior_facet_domains = max([-1] + [i.measure().domain_id() for i in form.interior_facet_integrals()]) + 1

    # Attach signature for convenience and reuse
    form_data.signature = form.signature()

    # Attach number of entries in element tensor
    dims = [create_element(v.element()).space_dimension() for v in form_data.arguments]
    dims_interior = [create_element(v.element()).space_dimension()*2 for v in form_data.arguments]
    form_data.num_entries = product(dims)
    form_data.num_entries_interior = product(dims_interior)

    # Attach number of facets
    form_data.num_facets = form_data.ffc_elements[0].cell().num_facets()

    end()
    return form, form_data

def compute_form_representations(form, form_data, options):
    "Compiler stage 2."

    begin("Compiler stage 2: Computing form representation(s)")
    representations = [Representation(form, form_data) for Representation in Representations]
    end()

    return representations

def optimize_form_representation(form, form_data):
    "Compiler stage 3."

    begin("Compiler stage 3: Optimizing form representation")
    info("Optimization of tensor contraction representation currently broken (to be fixed).")
    end()

def generate_form_code(form_data, representations, prefix, format, options):
    "Compiler stage 4."

    begin("Compiler stage 4: Generating code")

    # Generate common code like finite elements, dof map etc.
    common_code = generate_common_code(form_data, format)

    # Generate code for integrals
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
        else:
            forms.append(object)

    # Only compile element(s) when there are no forms
    if len(forms) > 0 and len(elements) > 0:
        elements = []

    return (forms, elements)

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

# FIXME: Old stuff below needs to be cleaned up
# FIXME: KBO: Is the above FIXME still valid? The function and class below are
# both used, and they look up to date and clean to me.
def create_ffc_coefficients(ufl_coefficients, ufl_coefficient_names):
    "Try to convert UFL functions to FFC Coefficients"

    class Coefficient:
        def __init__(self, element, name):
            self.element = element
            self.name = name
    ffc_coefficients = []
    for i, f in enumerate(ufl_coefficients):
        ffc_coefficients.append(Coefficient(create_element(f.element()), ufl_coefficient_names[i]))

    return ffc_coefficients

class ElementData:
    """This class holds meta data for a list of elements. It has the
    same attributes as FormData, but only the elements and their dof
    maps are available."""

    def __init__(self, elements):
        "Create element for list of elements"

        debug("Extracting element data...")

        self.signature                    = None
        self.rank                         = -1
        self.num_coefficients             = 0
        self.num_arguments                = len(elements)
        self.num_terms                    = 0
        self.num_cell_domains             = 0
        self.num_exterior_facet_domains   = 0
        self.num_interior_facet_domains   = 0
        self.elements                     = [create_element(element) for element in elements]
        self.dof_maps                     = [create_dof_map(element) for element in elements]
        self.ffc_elements                 = self.elements
        self.ffc_dof_maps                 = self.dof_maps
        self.coefficients                 = []

        debug("done")

# FFC format modules
from ufcformat import Format

