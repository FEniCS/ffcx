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
import copy
from itertools import chain

# UFL modules
from ufl.finiteelement import MixedElement, EnrichedElement
from ufl.algorithms import estimate_total_polynomial_degree
from ufl.algorithms import sort_elements
from ufl.algorithms import compute_form_data
from ufl.algorithms.analysis import extract_sub_elements

# FFC modules
from ffc.log import log, info, begin, end, warning, debug, error, ffc_assert, warning_blue
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
    return analyze_ufl_objects(forms, "form", parameters)


def analyze_elements(elements, parameters):
    return analyze_ufl_objects(elements, "element", parameters)


def analyze_coordinate_mappings(coordinate_elements, parameters):
    return analyze_ufl_objects(coordinate_elements, "coordinate_mapping", parameters)


def analyze_ufl_objects(ufl_objects, kind, parameters):
    """
    Analyze ufl object(s), either forms, elements, or coordinate mappings, returning:

       form_datas      - a tuple of form_data objects
       unique_elements - a tuple of unique elements across all forms
       element_numbers - a mapping to unique numbers for all elements

    """
    begin("Compiler stage 1: Analyzing %s(s)" % (kind,))

    form_datas = ()
    unique_elements = set()
    unique_coordinate_elements = set()

    if kind == "form":
        forms = ufl_objects

        # Analyze forms
        form_datas = tuple(_analyze_form(form, parameters)
                           for form in forms)

        # Extract unique elements accross all forms
        for form_data in form_datas:
            unique_elements.update(form_data.unique_sub_elements)

        # Extract coordinate elements across all forms
        for form_data in form_datas:
            unique_coordinate_elements.update(form_data.coordinate_elements)

    elif kind == "element":
        elements = ufl_objects

        # Extract unique (sub)elements
        unique_elements.update(extract_sub_elements(elements))

    elif kind == "coordinate_mapping":
        meshes = ufl_objects

        # Extract unique (sub)elements
        unique_coordinate_elements = [mesh.ufl_coordinate_element() for mesh in meshes]

    # Make sure coordinate elements and their subelements are included
    unique_elements.update(extract_sub_elements(unique_coordinate_elements))

    # Sort elements
    unique_elements = sort_elements(unique_elements)

    # Check for schemes for QuadratureElements
    for element in unique_elements:
        if element.family() == "Quadrature":
            qs = element.quadrature_scheme()
            if qs is None:
                error("Missing quad_scheme in quadrature element.")

    # Compute element numbers
    element_numbers = _compute_element_numbers(unique_elements)

    end()

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

    # Hack to override representation with environment variable
    forced_r = os.environ.get("FFC_FORCE_REPRESENTATION")
    if forced_r:
        warning("representation:    forced by $FFC_FORCE_REPRESENTATION to '%s'" % forced_r)

    # Compute form metadata
    if parameters["representation"] == "uflacs" or forced_r == "uflacs":
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


def _extract_common_quadrature_degree(integral_metadatas):
    # Check that quadrature degree is the same
    quadrature_degrees = [md["quadrature_degree"] for md in integral_metadatas]
    for d in quadrature_degrees:
        if not isinstance(d, int):
            error("Invalid non-integer quadrature degree %s" % (str(d),))
    qd = max(quadrature_degrees)
    if not all_equal(quadrature_degrees):
        # TODO: This may be loosened up without too much effort,
        # if the form compiler handles mixed integration degree,
        # something that most of the pipeline seems to be ready for.
        info("Quadrature degree must be equal within each sub domain, using degree %d." % qd)
    return qd


def _autoselect_quadrature_degree(integral_metadata, integral, form_data):
    # Automatic selection of quadrature degree
    qd = integral_metadata["quadrature_degree"]
    pd = integral_metadata["estimated_polynomial_degree"]

    # Special case: handling -1 as "auto" for quadrature_degree
    if qd == -1:
        qd = "auto"

    # TODO: Add other options here
    if qd == "auto":
        qd = pd
        info("quadrature_degree: auto --> %d" % qd)
    if isinstance(qd, int):
        if qd >= 0:
            info("quadrature_degree: %d" % qd)
        else:
            error("Illegal negative quadrature degree %s " % (qd,))
    else:
        error("Invalid quadrature_degree %s." % (qd,))

    tdim = integral.ufl_domain().topological_dimension()
    _check_quadrature_degree(qd, tdim)
    return qd


def _check_quadrature_degree(degree, top_dim):
    """Check that quadrature degree does not result in a unreasonable high
    number of integration points."""
    num_points = ((degree + 1 + 1) // 2)**top_dim
    if num_points >= 100:
        warning_blue("WARNING: The number of integration points for each cell will be: %d" % num_points)
        warning_blue("         Consider using the option 'quadrature_degree' to reduce the number of points")


def _extract_common_quadrature_rule(integral_metadatas):
    # Check that quadrature rule is the same
    # (To support mixed rules would be some work since num_points is
    #  used to identify quadrature rules in large parts of the pipeline)
    quadrature_rules = [md["quadrature_rule"] for md in integral_metadatas]
    if all_equal(quadrature_rules):
        qr = quadrature_rules[0]
    else:
        qr = "canonical"
        info("Quadrature rule must be equal within each sub domain, using %s rule." % qr)
    return qr


def _autoselect_quadrature_rule(integral_metadata, integral, form_data):
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
    # Return automatically determined quadrature rule
    return qr


def _determine_representation(integral_metadatas, ida, form_data, parameters):
    "Determine one unique representation considering all integrals together."

    # Hack because uflacs and quadrature/tensor cannot coincide in same form because of compute_form_data differences.
    r = parameters["representation"]
    if r == "uflacs":
        warning("representation:    ignoring metadata and using '%s' set by parameters" % r)
        return r

    # Hack to override representation with environment variable
    forced_r = os.environ.get("FFC_FORCE_REPRESENTATION")
    if forced_r:
        r = forced_r
        warning("representation:    forced by $FFC_FORCE_REPRESENTATION to '%s'" % r)
        return r

    # Check that representations are compatible
    # (Generating code with different representations within a
    # single tabulate_tensor is considered not worth the effort)
    representations = set()
    for md in integral_metadatas:
        if md["representation"] != "auto":
            representations.add(md["representation"])
    if len(representations) > 1:
        error("Integral representation must be equal within each sub domain or 'auto', got %s." % (str(list(set(representations))),))
    elif representations:
        r, = representations
    else:
        r = "auto"

    # If it's still auto, try to determine which representation is best for these integrals
    if r == "auto":
        rs = set()
        for integral in ida.integrals:
            rs.add(_auto_select_representation(integral,
                                               form_data.unique_sub_elements,
                                               form_data.function_replace_map))
        # If any failed to work with tensor, don't use tensor
        if "tensor" in rs and len(rs) > 1:
            rs.remove("tensor")
        # The end result must be unique
        if len(rs) != 1:
            error("Failed to auto-select representation, rs=%s." % (str(list(rs)),))
        r, = rs
        info("representation:    auto --> %s" % r)
    else:
        info("representation:    %s" % r)

    return r


def _attach_integral_metadata(form_data, parameters):
    "Attach integral metadata"
    # TODO: A nicer data flow would avoid modifying the form_data at all.

    # Recognized metadata keys
    metadata_keys = ("representation", "quadrature_degree", "quadrature_rule")
    metadata_parameters = {key: parameters[key] for key in metadata_keys if key in parameters}

    # Iterate over integral collections
    quad_schemes = []
    for ida in form_data.integral_data:
        # Iterate over integrals

        # Start with default values of integral metadata
        # (these will be either the FFC defaults, globally modified defaults,
        #  or overrides explicitly passed by the user to e.g. assemble())
        integral_metadatas = [copy.deepcopy(metadata_parameters)
                              for integral in ida.integrals]

        # Update with integral specific overrides
        for i, integral in enumerate(ida.integrals):
            integral_metadatas[i].update(integral.metadata() or {})

        # Determine representation, must be equal for all integrals on same subdomain
        r = _determine_representation(integral_metadatas, ida, form_data, parameters)
        for i, integral in enumerate(ida.integrals):
            integral_metadatas[i]["representation"] = r
        ida.metadata["representation"] = r

        # Determine automated updates to metadata values
        for i, integral in enumerate(ida.integrals):
            qr = _autoselect_quadrature_rule(integral_metadatas[i], integral, form_data)
            qd = _autoselect_quadrature_degree(integral_metadatas[i], integral, form_data)
            integral_metadatas[i]["quadrature_rule"] = qr
            integral_metadatas[i]["quadrature_degree"] = qd

        # Extract common metadata for integral collection
        qr = _extract_common_quadrature_rule(integral_metadatas)
        qd = _extract_common_quadrature_degree(integral_metadatas)
        ida.metadata["quadrature_rule"] = qr
        ida.metadata["quadrature_degree"] = qd

        # Add common num_cells (I believe this must be equal but I'm not that into this work)
        num_cells = set(md.get("num_cells") for md in integral_metadatas)
        if len(num_cells) != 1:
            error("Found integrals with different num_cells metadata on same subdomain: %s" % (str(list(num_cells)),))
        num_cells, = num_cells
        ida.metadata["num_cells"] = num_cells

        # Reconstruct integrals to avoid modifying the input integral,
        # which would affect the signature computation if the
        # integral was used again in the user program.
        # Modifying attributes of form_data.integral_data is less problematic
        # since it's lifetime is internal to the form compiler pipeline.
        for i, integral in enumerate(ida.integrals):
            ida.integrals[i] = integral.reconstruct(metadata=integral_metadatas[i])

        # Collect all quad schemes
        quad_schemes.extend([md["quadrature_rule"] for md in integral_metadatas])

    # Validate consistency of schemes for QuadratureElements
    # TODO: Can loosen up this a bit, only needs to be consistent
    # with the integrals that the elements are used in
    _validate_quadrature_schemes_of_elements(quad_schemes, form_data.unique_sub_elements)


def _validate_quadrature_schemes_of_elements(quad_schemes, elements):  # form_data):
    # Update scheme for QuadratureElements
    if quad_schemes and all_equal(quad_schemes):
        scheme = quad_schemes[0]
    else:
        scheme = "canonical"
        info("Quadrature rule must be equal within each sub domain, using %s rule." % scheme)
    for element in elements:
        if element.family() == "Quadrature":
            qs = element.quadrature_scheme()
            if qs != scheme:
                error("Quadrature element must have specified quadrature scheme (%s) equal to the integral (%s)." % (qs, scheme))


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
    if any(e.family() == "Quadrature" for e in sub_elements):
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
