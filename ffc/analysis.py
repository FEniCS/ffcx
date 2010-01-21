"""
Compiler stage 1: Analysis
--------------------------

This module implements the analysis/preprocessing of variational
forms, including automatic selection of elements, degrees and
form representation type.
"""

__author__ = "Anders Logg (logg@simula.no) and Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-02-05"
__copyright__ = "Copyright (C) 2007-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-21

# UFL modules
from ufl.common import istr
from ufl.finiteelement import MixedElement
from ufl.algorithms import preprocess, FormData
from ufl.algorithms import estimate_max_polynomial_degree
from ufl.algorithms import estimate_total_polynomial_degree
from ufl.algorithms import extract_unique_elements

# FFC modules
from ffc.log import log, info, begin, end
from ffc.quadratureelement import default_quadrature_degree

def analyze_forms(forms, object_names, options):
    """
    Analyze form(s), returning

       form_and_data   - a tuple of pairs (forms, form_data)
       unique_elements - a tuple of unique elements across all forms
       element_map     - a map from elements to unique element numbers
    """

    begin("Compiler stage 1: Analyzing form(s)")

    # Analyze forms
    form_and_data = [_analyze_form(form, object_names, options) for form in forms]

    # Extract unique elements
    unique_elements = []
    element_map = {}
    for (form, form_data) in form_and_data:
        for element in form_data.unique_sub_elements:
            if not element in element_map:
                element_map[element] = len(unique_elements)
                unique_elements.append(element)
    end()
    return form_and_data, unique_elements, element_map

def analyze_elements(elements):

    begin("Compiler stage 1: Analyzing form(s)")

    # Empty form and data
    form_and_data = []

    # FIXME: This looks unecessarily complex

    # Extract unique elements
    unique_elements = []
    element_map = {}
    for element in elements:
        # Check if element is present
        if not element in element_map:
            element_map[element] = len(unique_elements)
            unique_elements.append(element)
        # Check sub elements if any
        if isinstance(element, MixedElement):
            for sub_element in element.sub_elements():
                if not sub_element in element_map:
                    element_map[sub_element] = len(unique_elements)
                    unique_elements.append(sub_element)
    end()
    return form_and_data, unique_elements, element_map

def _analyze_form(form, object_names, options):
    "Analyze form, returning preprocessed form and form data."

    # Preprocess form
    if not form.is_preprocessed():
        form = preprocess(form)

    # Compute form data
    form_data = FormData(form, object_names=object_names)
    info(str(form_data))

    # Adjust cell and degree for elements when unspecified
    _adjust_elements(form_data)

    # Extract integral metadata
    form_data.metadata = _extract_metadata(form, options, form_data.elements)

    return form, form_data

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

def _extract_metadata(form, options, elements):
    "Check metadata for integral and return new integral with proper metadata."

    metadata = {}

    # Iterate over integrals
    for integral in form.integrals():

        # Set default values for metadata
        representation = options["representation"]
        quadrature_degree = options["quadrature_degree"]
        quadrature_rule = options["quadrature_rule"]

        if quadrature_rule is None:
            info("Quadrature rule: default")
        else:
            info("Quadrature rule: " + str(quadrature_rule))
        info("Quadrature degree: " + str(quadrature_degree))

        # Get metadata for integral (if any)
        integral_metadata = integral.measure().metadata() or {}
        for (key, value) in integral_metadata.iteritems():
            if key == "ffc_representation":
                representation = integral_metadata["ffc_representation"]
            elif key == "quadrature_degree":
                quadrature_degree = integral_metadata["quadrature_degree"]
            elif key == "quadrature_rule":
                quadrature_rule = integral_metadata["quadrature_rule"]
            else:
                warning("Unrecognized option '%s' for integral metadata." % key)

        # Check metadata
        valid_representations = ["tensor", "quadrature", "auto"]
        if not representation in valid_representations:
            error("Unrecognized form representation '%s', must be one of %s.",
                  representation, ", ".join("'%s'" % r for r in valid_representations))
        if quadrature_degree != "auto":
            try:
                quadrature_degree = int(quadrature_degree)
                if not quadrature_degree >= 0:
                    error("Illegal quadrature degree '%s' for integral, must be a nonnegative integer.",
                        str(quadrature_degree))
            except:
                error("Illegal quadrature degree '%s' for integral, must be a nonnegative integer or 'auto'.",
                    str(quadrature_degree))

        # Automatically select metadata if "auto" is selected
        if representation == "auto":
            representation = _auto_select_representation(integral)
        if quadrature_degree == "auto":
            quadrature_degree = _auto_select_quadrature_degree(integral, representation, elements)
        log(30, "Integral quadrature degree is %d." % quadrature_degree)

        # No quadrature rules have been implemented yet
        if quadrature_rule:
            warning("No quadrature rules have been implemented yet, using the default from FIAT.")

        # Set metadata for integral
        metadata[integral] = {"quadrature_degree": quadrature_degree,
                              "ffc_representation": representation,
                              "quadrature_rule":quadrature_rule}

    return metadata

def _auto_select_representation(integral):
    "Automatically select the best representation for integral."

    # FIXME: Implement this
    info("Automatic selection of representation not implemented, defaulting to quadrature.")
    return "tensor"

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
