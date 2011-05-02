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

# Modified by Marie E. Rognes, 2010

# Last changed: 2011-05-02

# UFL modules
from ufl.common import istr, tstr
from ufl.integral import Measure
from ufl.finiteelement import MixedElement, EnrichedElement
from ufl.algorithms import estimate_max_polynomial_degree
from ufl.algorithms import estimate_total_polynomial_degree
from ufl.algorithms import sort_elements
from ufl.algorithms import compute_form_arities

# FIXME: Import error when trying to import extract_sub_elements
# FIXME: from ufl.algorithmms.
from ufl.algorithms.analysis import extract_elements, extract_sub_elements

# FFC modules
from ffc.log import log, info, begin, end, warning, debug, error, ffc_assert
from ffc.utils import all_equal
from ffc.quadratureelement import default_quadrature_degree
from ffc.utils import all_equal
from ffc.tensor import estimate_cost

def analyze_forms(forms, object_names, parameters, common_cell=None):
    """
    Analyze form(s), returning

       form_datas      - a tuple of form_data objects
       unique_elements - a tuple of unique elements across all forms
       element_map     - a map from elements to unique element numbers
    """

    begin("Compiler stage 1: Analyzing form(s)")

    # Analyze forms
    form_datas = tuple(_analyze_form(form, object_names, parameters, common_cell) for form in forms)

    # Extract unique elements
    unique_elements = []
    for form_data in form_datas:
        for element in form_data.unique_sub_elements:
            if not element in unique_elements:
                unique_elements.append(element)

    # Sort elements
    unique_elements = sort_elements(unique_elements)

    # Build element map
    element_map = _build_element_map(unique_elements)

    end()

    return form_datas, unique_elements, element_map

def analyze_elements(elements, parameters):

    begin("Compiler stage 1: Analyzing form(s)")

    # Extract unique elements
    unique_elements = []
    element_map = {}
    for element in elements:
        # Get all (unique) nested elements.
        for e in _get_nested_elements(element):
            # Check if element is present
            if not e in element_map:
                element_map[e] = len(unique_elements)
                unique_elements.append(e)
    # Sort elements
    unique_elements = sort_elements(unique_elements)

    # Build element map
    element_map = _build_element_map(unique_elements)

    # Update scheme for QuadratureElements
    scheme = parameters["quadrature_rule"]
    if scheme == "auto":
        scheme = "default"
    for element in unique_elements:
        if element.family() == "Quadrature":
            element._quad_scheme = scheme
    end()

    return (), unique_elements, element_map

def _build_element_map(elements):
    "Build map from elements to element numbers."
    element_map = {}
    for (i, element) in enumerate(elements):
        element_map[element] = i
    return element_map

def _get_nested_elements(element):
    "Get unique nested elements (including self)."
    nested_elements = [element]
    for e in element.sub_elements():
        nested_elements += _get_nested_elements(e)
    return set(nested_elements)

def _analyze_form(form, object_names, parameters, common_cell=None):
    "Analyze form, returning preprocessed form."

    # Check that form is not empty
    ffc_assert(len(form.integrals()),
               "Form (%s) seems to be zero: cannot compile it." % str(form))

    # Extract element mapping for elements that should be replaced
    element_mapping, common_cell = _extract_element_mapping(form, common_cell)

    # Compute form metadata
    form_data = form.compute_form_data(object_names, common_cell, element_mapping)
    info("")
    info(str(form_data))

    # Extract preprocessed form
    preprocessed_form = form_data.preprocessed_form

    # Check that all terms in form have same arity
    ffc_assert(len(compute_form_arities(preprocessed_form)) == 1,
               "All terms in form must have same rank.")

    # Adjust cell and degree for elements when unspecified
    # FIXME: FiniteElementBase.set_foo() will not be supported from UFL 1.0.
    _adjust_elements(form_data)

    # FIXME: Does this function also modify the form directly and should
    # FIXME: should therefore be removed (like the UFL set_foo functions)
    # Extract integral metadata
    _extract_metadata(form_data, parameters)

    return form_data

def _extract_element_mapping(form, common_cell):
    """Extract mapping for elements that should be replaced. This is
    used to automatically set elements in cases where the element is
    not completely specified, which is typically the case for DOLFIN
    Expressions."""

    print "Extracting element mapping"

    # Extract all elements
    elements = extract_elements(form)
    elements = extract_sub_elements(elements)

    # Store cell
    if common_cell is None:
        cells = [e.cell() for e in elements]
        cells = [c for c in cells if c.domain() is not None]
        if len(cells) == 0:
            error("Unable to extract common element; missing cell definition in form.")
        common_cell = cells[0]

    # FIXME: Check if these are needed
    #    elif form._integrals:
    #        # Special case to allow functionals only depending on geometric variables, with no elements
    #        form_data.cell = form._integrals[0].integrand().cell()
    #    else:
    #        # Special case to allow integral of constants to pass through without crashing
    #        form_data.cell = None
    #        warning("Form is empty, no elements or integrals, cell is undefined.")

    # Extract common degree
    common_degree = max([e.degree() for e in elements])
    if common_degree is None:
        common_degree = default_quadrature_degree

    # Degree must be at least 1 (to work with Lagrange elements)
    common_degree = max(1, common_degree)

    # Reconstruct elements
    element_mapping = {}
    for element in elements:

        # Adjust element family (currently unchanged)
        family = element.family()

        # Adjust cell
        cell = element.cell()
        if cell.domain() is None:
            info("Adjusting element cell from %s to %s." % (istr(cell), str(common_cell)))
            cell = common_cell

        # Adjust degree
        degree = element.degree()
        if degree is None:
            info("Adjusting element degree from %s to %d" % (istr(degree), common_degree))
            degree = common_degree

        # Reconstruct element and set mapping
        #new_element = element.reconstruct(family=family, cell=cell, degree=degree)
        #element_mapping[element] = new_element

    return element_mapping, common_cell

def _adjust_elements(form_data):
    "Adjust cell and degree for elements when unspecified."

    # Note the importance of consider form_data.sub_elements here
    # instead of form_data.unique_sub_elements. This is because
    # elements considered equal (same output from __repr__) will not
    # be repeated in unique_sub_elements but all elements need to be
    # adjusted.

    # Extract common cell
    common_cell = form_data.cell
    if common_cell.domain() is None:
        error("Missing cell definition in form.")

    # Extract common degree
    common_degree = max([element.degree() for element in form_data.sub_elements])
    if common_degree is None:
        common_degree = default_quadrature_degree

    # Degree must be at least 1 (to work with Lagrange elements)
    common_degree = max(1, common_degree)

    # Set cell and degree if missing
    # FIXME: FiniteElementBase.set_foo() will not be supported from UFL 1.0.
    for element in form_data.sub_elements:

        # Check if cell and degree need to be adjusted
        cell = element.cell()
        degree = element.degree()
        if degree is None:
            info("Adjusting element degree from %s to %d" % (istr(degree), common_degree))
            element.set_degree(common_degree)
        if cell.domain() is None:
            info("Adjusting element cell from %s to %s." % (istr(cell), str(common_cell)))
            element.set_cell(common_cell)

def _extract_metadata(form_data, parameters):
    "Attach and group meta data for each subdomain integral collection."

    # Recognized metadata keys
    metadata_keys = ("representation", "quadrature_degree", "quadrature_rule")

    quad_schemes = []
    # Iterate over integral collections
    for (domain_type, domain_id, integrals, metadata) in form_data.integral_data:

        # Iterate over integrals
        integral_metadatas = []
        for integral in integrals:

            # Get metadata for integral
            integral_metadata = integral.measure().metadata() or {}
            for key in metadata_keys:
                if not key in integral_metadata:
                    integral_metadata[key] = parameters[key]

            # Check metadata
            r  = integral_metadata["representation"]
            qd = integral_metadata["quadrature_degree"]
            qr = integral_metadata["quadrature_rule"]
            if not r in ("quadrature", "tensor", "auto"):
                info("Valid choices are 'tensor', 'quadrature' or 'auto'.")
                error("Illegal choice of representation for integral: " + str(r))
            if not qd  == "auto":
                qd = int(qd)
                if not qd >= 0:
                    info("Valid choices are nonnegative integers or 'auto'.")
                    error("Illegal quadrature degree for integral: " + str(qd))
            if not qr in ("default", "canonical", "auto"):
                info("Valid choices are 'default', 'canonical' or 'auto'.")
                error("Illegal choice of quadrature rule for integral: " + str(qr))

            # Automatic selection of representation
            if r == "auto":
                r = _auto_select_representation(integral, form_data.unique_sub_elements)
                info("representation:    auto --> %s" % r)
                integral_metadata["representation"] = r
            else:
                info("representation:    %s" % r)

            # Automatic selection of quadrature degree
            if qd == "auto":
                qd = _auto_select_quadrature_degree(integral, r, form_data.unique_sub_elements)
                info("quadrature_degree: auto --> %d" % qd)
                integral_metadata["quadrature_degree"] = qd
            else:
                info("quadrature_degree: %d" % qd)

            # Automatic selection of quadrature rule
            if qr == "auto":
                # Just use default for now.
                qr = "default"
                info("quadrature_rule:   auto --> %s" % qr)
                integral_metadata["quadrature_rule"] = qr
            else:
                info("quadrature_rule:   %s" % qr)
            quad_schemes.append(qr)

            # Append to list of metadata
            integral_metadatas.append(integral_metadata)

        # Extract common metadata for integral collection
        if len(integrals) == 1:
            metadata.update(integral_metadatas[0])
        else:

            # Check that representation is the same
            # FIXME: Why must the representation within a sub domain be the same?
            representations = [md["representation"] for md in integral_metadatas]
            if not all_equal(representations):
                r = "quadrature"
                info("Integral representation must be equal within each sub domain, using %s representation." % r)
            else:
                r = representations[0]

            # Check that quadrature degree is the same
            # FIXME: Why must the degree within a sub domain be the same?
            quadrature_degrees = [md["quadrature_degree"] for md in integral_metadatas]
            if not all_equal(quadrature_degrees):
                qd = max(quadrature_degrees)
                info("Quadrature degree must be equal within each sub domain, using degree %d." % qd)
            else:
                qd = quadrature_degrees[0]

            # Check that quadrature rule is the same
            # FIXME: Why must the rule within a sub domain be the same?
            quadrature_rules = [md["quadrature_rule"] for md in integral_metadatas]
            if not all_equal(quadrature_rules):
                qr = "canonical"
                info("Quadrature rule must be equal within each sub domain, using %s rule." % qr)
            else:
                qr = quadrature_rules[0]

            # Update common metadata
            metadata["representation"] = r
            metadata["quadrature_degree"] = qd
            metadata["quadrature_rule"] = qr

    # Update scheme for QuadratureElements
    if not all_equal(quad_schemes):
        scheme = "canonical"
        info("Quadrature rule must be equal within each sub domain, using %s rule." % qr)
    else:
        scheme = quad_schemes[0]
    for element in form_data.sub_elements:
        if element.family() == "Quadrature":
            element._quad_scheme = scheme

def _get_sub_elements(element):
    "Get sub elements."
    sub_elements = [element]
    if isinstance(element, MixedElement):
        for e in element.sub_elements():
            sub_elements += _get_sub_elements(e)
    elif isinstance(element, EnrichedElement):
        for e in element._elements:
            sub_elements += _get_sub_elements(e)
    return sub_elements

def _auto_select_representation(integral, elements):
    """
    Automatically select a suitable representation for integral.
    Note that the selection is made for each integral, not for
    each term. This means that terms which are grouped by UFL
    into the same integral (if their measures are equal) will
    necessarily get the same representation.
    """

    # Get ALL sub elements, needed to check for restrictions of EnrichedElements.
    sub_elements = []
    for e in elements:
        sub_elements += _get_sub_elements(e)

    # Use quadrature representation if we have a quadrature element
    if len([e for e in sub_elements if e.family() == "Quadrature"]):
        return "quadrature"

    # Use quadrature representation if any elements are restricted to
    # UFL.Measure. This is used when integrals are computed over discontinuities.
    if len([e for e in sub_elements if isinstance(e.domain_restriction(), Measure)]):
        return "quadrature"

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
            debug("Found quadrature element(s) with the following degree(s): " + str(quadrature_degrees))
            ffc_assert(min(quadrature_degrees) == max(quadrature_degrees), \
                       "All QuadratureElements in an integrand must have the same degree: %s" \
                       % str(quadrature_degrees))
            debug("Selecting quadrature degree based on quadrature element: " + str(quadrature_degrees[0]))
            return quadrature_degrees[0]

    # Otherwise estimate total degree of integrand
    q = estimate_total_polynomial_degree(integral, default_quadrature_degree)
    debug("Selecting quadrature degree based on total polynomial degree of integrand: " + str(q))

    return q
