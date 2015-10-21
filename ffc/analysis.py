"""
Compiler stage 1: Analysis
--------------------------

This module implements the analysis/preprocessing of variational
forms, including automatic selection of elements, degrees and
form representation type.
"""

# Copyright (C) 2007-201r Anders Logg and Kristian B. Oelgaard
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Marie E. Rognes, 2010
# Modified by Martin Alnaes, 2013-2014

import os
from itertools import chain

# UFL modules
from ufl.finiteelement import MixedElement, EnrichedElement
from ufl.algorithms import estimate_total_polynomial_degree
from ufl.algorithms import sort_elements
from ufl.algorithms import compute_form_data
from ufl.algorithms.analysis import extract_sub_elements

# FFC modules
from ffc.log import log, info, begin, end, warning, debug, error, ffc_assert, warning_blue
from ffc.utils import all_equal
from ffc.quadratureelement import default_quadrature_degree
from ffc.utils import all_equal
from ffc.tensor import estimate_cost

def analyze_forms(forms, parameters):
    """
    Analyze form(s), returning

       form_datas      - a tuple of form_data objects
       unique_elements - a tuple of unique elements across all forms
       element_numbers - a mapping to unique numbers for all elements
    """

    begin("Compiler stage 1: Analyzing form(s)")

    # Analyze forms
    form_datas = tuple(_analyze_form(form,
                                     parameters) for form in forms)

    # Extract unique elements accross all forms
    unique_elements = set()
    for form_data in form_datas:
        unique_elements.update(form_data.unique_sub_elements)

    # Sort elements
    unique_elements = sort_elements(unique_elements)

    # Compute element numbers
    element_numbers = _compute_element_numbers(unique_elements)

    # Extract coordinate elements
    unique_coordinate_elements = sorted(set(chain(*[form_data.coordinate_elements for form_data in form_datas])))

    end()

    return form_datas, unique_elements, element_numbers, unique_coordinate_elements

def analyze_elements(elements, parameters):

    begin("Compiler stage 1: Analyzing form(s)")

    # Extract unique (sub)elements
    unique_elements = set(extract_sub_elements(elements))

    # Sort elements
    unique_elements = sort_elements(unique_elements)

    # Build element map
    element_numbers = _compute_element_numbers(unique_elements)

    # Update scheme for QuadratureElements
    scheme = parameters["quadrature_rule"]
    if scheme == "auto":
        scheme = "default"
    for element in unique_elements:
        if element.family() == "Quadrature":
            element._quad_scheme = scheme

    end()

    form_datas = ()
    unique_coordinate_elements = ()
    return form_datas, unique_elements, element_numbers, unique_coordinate_elements

def _compute_element_numbers(elements):
    "Build map from elements to element numbers."
    element_numbers = {}
    for (i, element) in enumerate(elements):
        element_numbers[element] = i
    return element_numbers

def _analyze_form(form, parameters):
    "Analyze form, returning form data."

    # Check that form is not empty
    ffc_assert(not form.empty(),
               "Form (%s) seems to be zero: cannot compile it." % str(form))

    # Compute form metadata
    if parameters["representation"] == "uflacs":
        # Temporary workaround to let uflacs have a different preprocessing pipeline
        # than the legacy representations quadrature and tensor. This approach imposes
        # a limitation that e.g. uflacs and tensor representation cannot be mixed in the same form.
        from ufl.classes import Jacobian
        form_data = compute_form_data(form,
                                      do_apply_function_pullbacks=True,
                                      do_apply_integral_scaling=True,
                                      do_apply_geometry_lowering=True,
                                      preserve_geometry_types=(Jacobian,),
                                      do_apply_restrictions=True,
                                      )
    else:
        form_data = compute_form_data(form)

    info("")
    info(str(form_data))

    # Attach integral meta data
    _attach_integral_metadata(form_data, parameters)

    return form_data

# FIXME: Refactor this code. The data flow here is really something special. It's also buggy with side effects into user code by modifying dicts in-place.
def _attach_integral_metadata(form_data, parameters):
    "Attach integral metadata"

    # Recognized metadata keys
    metadata_keys = ("representation", "quadrature_degree", "quadrature_rule")

    # Iterate over integral collections
    quad_schemes = []
    for ida in form_data.integral_data:
        # TODO: Is it possible to detach this from IntegralData? It's a bit strange from the ufl side.
        common_metadata = ida.metadata

        # Iterate over integrals
        integral_metadatas = []
        for integral in ida.integrals:

            # Fill in integral metadata with default values
            # NB! This modifies the metadata of the input integral data!
            integral_metadata = integral.metadata() or {}
            for key in metadata_keys:
                if key not in integral_metadata:
                    integral_metadata[key] = parameters[key]

            # Automatic selection of representation
            r = integral_metadata["representation"]

            # Hack to override representation with environment variable
            forced_r = os.environ.get("FFC_FORCE_REPRESENTATION")
            if forced_r:
                r = forced_r
                warning("representation:    forced by $FFC_FORCE_REPRESENTATION to '%s'" % r)
            elif r == "auto":
                r = _auto_select_representation(integral,
                                                form_data.unique_sub_elements,
                                                form_data.function_replace_map)
                info("representation:    auto --> %s" % r)
            elif r in ("quadrature", "tensor", "uflacs"):
                info("representation:    %s" % r)
            else:
                info("Valid choices are 'tensor', 'quadrature', 'uflacs', or 'auto'.")
                error("Illegal choice of representation for integral: " + str(r))
            integral_metadata["representation"] = r

            # Automatic selection of quadrature degree
            qd = integral_metadata["quadrature_degree"]
            # Special case: handling -1 as "auto" for quadrature_degree
            if qd in ("auto", -1):
                qd = _auto_select_quadrature_degree(integral.integrand(),
                                                    r,
                                                    form_data.unique_sub_elements,
                                                    form_data.element_replace_map)
                info("quadrature_degree: auto --> %d" % qd)
            else:
                qd = int(qd)
                info("quadrature_degree: %d" % qd)
            # Validate degree
            if not qd >= 0:
                info("Valid choices are nonnegative integers or 'auto'.")
                error("Illegal quadrature degree for integral: " + str(qd))
            tdim = integral.ufl_domain().topological_dimension()
            _check_quadrature_degree(qd, tdim)

            integral_metadata["quadrature_degree"] = qd
            assert isinstance(qd, int)

            # Automatic selection of quadrature rule
            qr = integral_metadata["quadrature_rule"]
            if qr == "auto":
                # Just use default for now.
                qr = "default"
                info("quadrature_rule:   auto --> %s" % qr)
            elif qr in ("default", "canonical", "vertex"):
                info("quadrature_rule:   %s" % qr)
            else:
                info("Valid choices are 'default', 'canonical', 'vertex', and 'auto'.")
                error("Illegal choice of quadrature rule for integral: " + str(qr))
            integral_metadata["quadrature_rule"] = qr
            quad_schemes.append(qr)

            # Append to list of metadata
            integral_metadatas.append(integral_metadata)

        # Extract common metadata for integral collection
        if len(ida.integrals) == 1:
            common_metadata.update(integral_metadatas[0])
        else:
            # Check that representation is the same
            # (Generating code with different representations within a
            # single tabulate_tensor is considered not worth the effort)
            representations = [md["representation"] for md in integral_metadatas]
            if all_equal(representations):
                r = representations[0]
            else:
                r = "quadrature"
                info("Integral representation must be equal within each sub domain, using %s representation." % r)

            # Check that quadrature degree is the same
            # FIXME: Why must the degree within a sub domain be the same?
            #        This makes no sense considering that num_points is
            #        used as a key all over in quadrature representation...
            quadrature_degrees = [md["quadrature_degree"] for md in integral_metadatas]
            if all_equal(quadrature_degrees):
                qd = quadrature_degrees[0]
            else:
                qd = max(quadrature_degrees)
                info("Quadrature degree must be equal within each sub domain, using degree %d." % qd)
            assert isinstance(qd, int)

            # Check that quadrature rule is the same
            # FIXME: Why must the rule within a sub domain be the same?
            #        To support this would be more work since num_points is used
            #        to identify quadrature rules in the quadrature representation.
            quadrature_rules = [md["quadrature_rule"] for md in integral_metadatas]
            if all_equal(quadrature_rules):
                qr = quadrature_rules[0]
            else:
                qr = "canonical"
                info("Quadrature rule must be equal within each sub domain, using %s rule." % qr)

            # Update common metadata
            assert isinstance(qd, int)
            common_metadata["representation"] = r
            common_metadata["quadrature_degree"] = qd
            common_metadata["quadrature_rule"] = qr

    # Update scheme for QuadratureElements
    if quad_schemes and all_equal(quad_schemes):
        scheme = quad_schemes[0]
    else:
        scheme = "canonical"
        info("Quadrature rule must be equal within each sub domain, using %s rule." % scheme)

    # FIXME: This modifies the elements depending on the form compiler parameters,
    #        this is a serious breach of the immutability of ufl objects, since the
    #        element quad scheme is part of the signature and hash of the element...
    for element in form_data.unique_sub_elements:
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

def _auto_select_representation(integral, elements, function_replace_map):
    """
    Automatically select a suitable representation for integral.
    Note that the selection is made for each integral, not for
    each term. This means that terms which are grouped by UFL
    into the same integral (if their measures are equal) will
    necessarily get the same representation.
    """

    # Skip unsupported integration domain types
    if integral.integral_type() == "vertex":
        return "quadrature"

    # Get ALL sub elements, needed to check for restrictions of EnrichedElements.
    sub_elements = []
    for e in elements:
        sub_elements += _get_sub_elements(e)

    # Use quadrature representation if we have a quadrature element
    if len([e for e in sub_elements if e.family() == "Quadrature"]):
        return "quadrature"

    # Estimate cost of tensor representation
    tensor_cost = estimate_cost(integral, function_replace_map)
    debug("Estimated cost of tensor representation: " + str(tensor_cost))

    # Use quadrature if tensor representation is not possible
    if tensor_cost == -1:
        return "quadrature"

    # Otherwise, select quadrature when cost is high
    if tensor_cost <= 3:
        return "tensor"
    else:
        return "quadrature"

def _auto_select_quadrature_degree(integrand, representation, elements, element_replace_map):
    "Automatically select a suitable quadrature degree for integrand."
    # TODO: Move this to form preprocessing, as part of integral_data?

    # Use quadrature element degree if any is found
    quadrature_degrees = [e.degree() for e in elements if e.family() == "Quadrature"]
    if quadrature_degrees:
        debug("Found quadrature element(s) with the following degree(s): " + str(quadrature_degrees))
        ffc_assert(min(quadrature_degrees) == max(quadrature_degrees), \
                   "All QuadratureElements in an integrand must have the same degree: %s" \
                   % str(quadrature_degrees))
        debug("Selecting quadrature degree based on quadrature element: " + str(quadrature_degrees[0]))
        ffc_assert(representation != "tensor", "Tensor representation does not support quadrature elements.")
        return quadrature_degrees[0]

    # Otherwise estimate total degree of integrand
    q = estimate_total_polynomial_degree(integrand, default_quadrature_degree, element_replace_map)
    debug("Selecting quadrature degree based on total polynomial degree of integrand: " + str(q))

    return q

def _check_quadrature_degree(degree, top_dim):
    """Check that quadrature degree does not result in a unreasonable high
    number of integration points."""
    num_points = ((degree + 1 + 1) // 2)**top_dim
    if num_points >= 100:
        warning_blue("WARNING: The number of integration points for each cell will be: %d" % num_points)
        warning_blue("         Consider using the option 'quadrature_degree' to reduce the number of points")
