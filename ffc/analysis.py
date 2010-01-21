"""
Compiler stage 1: Analysis
--------------------------

This module implements the analysis/preprocessing of variational
forms, including automatic selection of elements, degrees and
form representation type.
"""

__author__ = "Anders Logg (logg@simula.no) and Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2007-02-05"
__copyright__ = "Copyright (C) 2007-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-21

# UFL modules
from ufl.common import istr, tstr
from ufl.finiteelement import MixedElement
from ufl.algorithms import preprocess, FormData
from ufl.algorithms import estimate_max_polynomial_degree
from ufl.algorithms import estimate_total_polynomial_degree
from ufl.algorithms import extract_unique_elements

# FFC modules
from ffc.log import log, info, begin, end, warning, debug
from ffc.quadratureelement import default_quadrature_degree
from ffc.tensor import estimate_cost

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
    "Adjust cell and degree for elements when unspecified."

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
            info("Adjusting element degree from %s to %d" % (istr(degree), common_degree))
            element.set_degree(common_degree)
        if cell.domain() is None:
            info("Adjusting element cell from %s to %s." % (istr(cell), str(common_cell)))
            element.set_cell(common_cell)

def _extract_metadata(form, options, elements):
    "Check metadata for integral and return new integral with proper metadata."

    metadata = {}

    # Recognized metadata keys
    metadata_keys = ("representation", "quadrature_degree")

    # Iterate over integrals
    for integral in form.integrals():

        # Get metadata for integral
        integral_metadata = integral.measure().metadata() or {}
        for key in metadata_keys:
            if not key in integral_metadata:
                integral_metadata[key] = options[key]

        # Check metadata
        r = integral_metadata["representation"]
        q = integral_metadata["quadrature_degree"]
        if not r in ("quadrature", "tensor", "auto"):
            info("Valid choices are 'tensor', 'quadrature' or 'auto'.")
            error("Illegal choice of representation for integral: " + str(r))
        if not ((isinstance(q, int) and q < 0) or q == "auto"):
            info("Valid choices are nonnegative integers or 'auto'.")
            error("Illegal quadrature degree for integral: " + str(q))

        # Automatic selection of representation
        if r == "auto":
            r = _auto_select_representation(integral)
            info("representation:    auto --> %s" % r)
            integral_metadata["representation"] = r
        else:
            info("representation:    %s" % s)

        # Automatic selection of quadrature degree
        if q == "auto":
            q = _auto_select_quadrature_degree(integral, r, elements)
            info("quadrature_degree: auto --> %d" % q)
            integral_metadata["quadrature_degree"] = q
        else:
            info("quadrature_degree: %d" % q)

        info("")

        # Quadrature rule selection not supported
        integral_metadata["quadrature_rule"] = "auto"

        # Set metadata for integral
        metadata[integral] = integral_metadata

    return metadata

def _auto_select_representation(integral):
    "Automatically select a suitable representation for integral."

    # Estimate cost of tensor representation
    tensor_cost = estimate_cost(integral)
    debug("Estimated cost of tensor representation: " + str(tensor_cost))

    # Use quadrature if tensor representation is not possible
    if tensor_cost == -1:
        return "quadrature"

    # Otherwise, select quadrature when cost is high
    if tensor_cost <= 3:
        return "tensor"
    else:
        return "quadrature"

def _auto_select_quadrature_degree(integral, representation, elements):
    "Automatically select a suitable quadrature degree for integral."

    # Use maximum quadrature element degree if any for quadrature representation
    if representation == "quadrature":
        quadrature_degrees = [e.degree() for e in elements if e.family() == "Quadrature"]
        if quadrature_degrees:
            ffc_assert(min(quadrature_degrees) == max(quadrature_degrees), \
                       "All QuadratureElements in an integrand must have the same degree: %s" \
                       % str(quadrature_degrees))
            return quadrature_degrees[0]

    # Otherwise estimate total degree of integrand
    return estimate_total_polynomial_degree(integral, default_quadrature_degree)
