# -*- coding: utf-8 -*-

# Copyright (C) 2007-2017 Anders Logg, Martin Alnaes, Kristian B. Oelgaard,
# and others
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

"""
Compiler stage 1: Analysis
--------------------------

This module implements the analysis/preprocessing of variational
forms, including automatic selection of elements, degrees and
form representation type.
"""

import numpy
import os
import copy
from itertools import chain

# UFL modules
from ufl.classes import Form, CellVolume, FacetArea
from ufl.integral import Integral
from ufl.finiteelement import MixedElement, EnrichedElement, VectorElement
from ufl.algorithms import sort_elements
from ufl.algorithms import compute_form_data
from ufl.algorithms.analysis import extract_sub_elements
from ufl import custom_integral_types

# FFC modules
from ffc.log import info, begin, end, warning, debug, error, ffc_assert, warning_blue
from ffc.utils import all_equal

# Default precision for formatting floats
default_precision = numpy.finfo("double").precision + 1  # == 16

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
    #unique_coordinate_elements = sort_elements(unique_coordinate_elements)
    unique_coordinate_elements = sorted(unique_coordinate_elements, key=lambda x: repr(x))

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
    if form.empty():
        error("Form (%s) seems to be zero: cannot compile it." % str(form))

    # Hack to override representation with environment variable
    forced_r = os.environ.get("FFC_FORCE_REPRESENTATION")
    if forced_r:
        warning("representation:    forced by $FFC_FORCE_REPRESENTATION to '%s'" % forced_r)
        r = forced_r
    else:
        # Check representation parameters to figure out how to
        # preprocess
        r = _extract_representation_family(form, parameters)
    debug("Preprocessing form using '%s' representation family." % r)

    # Compute form metadata
    if r == "uflacs":
        # Temporary workaround to let uflacs have a different
        # preprocessing pipeline than the legacy quadrature
        # representation. This approach imposes a limitation that,
        # e.g. uflacs and qudrature, representations cannot be mixed
        # in the same form.
        from ufl.classes import Jacobian
        form_data = compute_form_data(form,
                                      do_apply_function_pullbacks=True,
                                      do_apply_integral_scaling=True,
                                      do_apply_geometry_lowering=True,
                                      preserve_geometry_types=(Jacobian,),
                                      do_apply_restrictions=True)
    elif r == "tsfc":
        try:
            # TSFC provides compute_form_data wrapper using correct
            # kwargs
            from tsfc.ufl_utils import compute_form_data as tsfc_compute_form_data
        except ImportError:
            error("Failed to import tsfc.ufl_utils.compute_form_data when asked "
                  "for tsfc representation.")
        form_data = tsfc_compute_form_data(form)
    elif r == "quadrature":
        # quadrature representation
        form_data = compute_form_data(form)
    else:
        error("Unexpected representation family '%s' for form preprocessing." % r)

    info("")
    info(str(form_data))

    # Attach integral meta data
    _attach_integral_metadata(form_data, r, parameters)
    _validate_representation_choice(form_data, r)

    return form_data


def _extract_representation_family(form, parameters):
    """Return 'uflacs', 'tsfc' or 'quadrature', or raise error. This takes
    care of (a) compatibility between representations due to
    differences in preprocessing, (b) choosing uflacs for higher-order
    geometries.

    NOTE: Final representation is picked later by
    ``_determine_representation``.

    """

    # Fetch all representation choice requests from metadata
    representations = set(integral.metadata().get("representation", "auto")
                          for integral in form.integrals())

    # If auto is present, add parameters value (which may still be
    # auto) and then remove auto so there's no auto in the set
    if "auto" in representations:
        representations.add(parameters["representation"])
        representations.remove("auto")

    # Sanity check
    ffc_assert(len(representations.intersection(('auto', None))) == 0,
               "Unexpected representation family candidates '%s'." % representations)

    # No representations requested, find compatible representations
    compatible = _find_compatible_representations(form.integrals(), [])

    if len(representations) == 1:
        r = representations.pop()
        if r not in compatible:
            error("Representation family %s is not compatible with this form (try one of %s)" % (r, sorted(compatible)))
        return r
    elif len(representations) == 0:
        if len(compatible) == 1:
            # If only one compatible, use it
            return compatible.pop()
        else:
            # Default to uflacs
            # NOTE: Need to pick the same default as in _auto_select_representation
            return "uflacs"
    else:
        # Don't tolerate user requests for mixing old and new
        # representation families in same form due to restrictions in
        # preprocessing
        assert len(representations) > 1
        error("Cannot mix quadrature, uflacs, or tsfc "
              "representation in single form.")


def _validate_representation_choice(form_data,
                                    preprocessing_representation_family):
    """Check that effective representations

    * do not mix quadrature, uflacs and tsfc,
    * implement higher-order geometry,
    * match employed preprocessing strategy.

    This function is final check that everything is compatible due
    to the mess in this file. Better safe than sorry...

    """

    # Fetch all representations from integral metadata (this should by
    # now be populated with parameters or defaults instead of auto)
    representations = set()
    for ida in form_data.integral_data:
        representations.add(ida.metadata["representation"])
        for integral in ida.integrals:
            representations.add(integral.metadata()["representation"])

    # No integrals
    if len(representations) == 0:
        return

    # Require unique family; allow quadrature only with affine meshes
    if len(representations) != 1:
        error("Failed to extract unique representation family. "
              "Got '%s'." % representations)

    if _has_higher_order_geometry(form_data.preprocessed_form):
        ffc_assert('quadrature' not in representations,
            "Did not expect quadrature representation for higher-order geometry.")

    # Check preprocessing strategy
    ffc_assert(preprocessing_representation_family in representations,
               "Form has been preprocessed using '%s' representaion family, "
               "while '%s' representations have been set for integrals."
               % (preprocessing_representation_family, representations))


def _has_custom_integrals(o):
    if isinstance(o, Integral):
        return o.integral_type() in custom_integral_types
    elif isinstance(o, Form):
        return any(_has_custom_integrals(itg) for itg in o.integrals())
    elif isinstance(o, (list, tuple)):
        return any(_has_custom_integrals(itg) for itg in o)
    else:
        raise NotImplementedError


def _has_higher_order_geometry(o):
    if isinstance(o, Integral):
        P1 = VectorElement("P", o.ufl_domain().ufl_cell(), 1)
        return o.ufl_domain().ufl_coordinate_element() != P1
    elif isinstance(o, Form):
        P1 = VectorElement("P", o.ufl_cell(), 1)
        return any(d.ufl_coordinate_element() != P1 for d in o.ufl_domains())
    elif isinstance(o, (list, tuple)):
        return any(_has_higher_order_geometry(itg) for itg in o)
    else:
        raise NotImplementedError


def _extract_common_quadrature_degree(integral_metadatas):
    # Check that quadrature degree is the same
    quadrature_degrees = [md["quadrature_degree"] for md in integral_metadatas]
    for d in quadrature_degrees:
        if not isinstance(d, int):
            error("Invalid non-integer quadrature degree %s" % (str(d),))
    qd = max(quadrature_degrees)
    if not all_equal(quadrature_degrees):
        # FIXME: Shouldn't we raise here?
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
    if qd in [-1, None]:
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
    # Check that quadrature rule is the same (To support mixed rules
    # would be some work since num_points is used to identify
    # quadrature rules in large parts of the pipeline)
    quadrature_rules = [md["quadrature_rule"] for md in integral_metadatas]
    if all_equal(quadrature_rules):
        qr = quadrature_rules[0]
    else:
        qr = "canonical"
        # FIXME: Shouldn't we raise here?
        info("Quadrature rule must be equal within each sub domain, using %s rule." % qr)
    return qr


def _autoselect_quadrature_rule(integral_metadata, integral, form_data):
    # Automatic selection of quadrature rule
    qr = integral_metadata["quadrature_rule"]
    if qr in ["auto", None]:
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


def _determine_representation(integral_metadatas, ida, form_data, form_r_family,
                              parameters):
    "Determine one unique representation considering all integrals together."

    # Extract unique representation among these single-domain
    # integrals (Generating code with different representations within
    # a single tabulate_tensor is considered not worth the effort)
    representations  = set(md["representation"] for md in integral_metadatas
                           if md["representation"] != "auto")
    optimize_values  = set(md["optimize"] for md in integral_metadatas)
    precision_values = set(md["precision"] for md in integral_metadatas)

    if len(representations) > 1:
        error("Integral representation must be equal within each sub domain or 'auto', got %s." % (str(sorted(str(v) for v in representations)),))
    if len(optimize_values) > 1:
        error("Integral 'optimize' metadata must be equal within each sub domain or not set, got %s." % (str(sorted(str(v) for v in optimize_values)),))
    if len(precision_values) > 1:
        error("Integral 'precision' metadata must be equal within each sub domain or not set, got %s." % (str(sorted(str(v) for v in precision_values)),))

    # The one and only non-auto representation found, or get from parameters
    r, = representations  or (parameters["representation"],)
    o, = optimize_values  or (parameters["optimize"],)
    # FIXME: Default param value is zero which is not interpreted well by tsfc!
    p, = precision_values or (parameters["precision"],)

    # If it's still auto, try to determine which representation is
    # best for these integrals
    if r == "auto":
        # Find representations compatible with these integrals
        compatible = _find_compatible_representations(ida.integrals,
                                                      form_data.unique_sub_elements)
        # Pick the one compatible or default to uflacs
        if len(compatible) == 0:
            error("Found no representation capable of compiling this form.")
        elif len(compatible) == 1:
            r, = compatible
        else:
            # NOTE: Need to pick the same default as in
            # _extract_representation_family
            if form_r_family == "uflacs":
                r = "uflacs"
            elif form_r_family == "tsfc":
                r = "tsfc"
            elif form_r_family == "quadrature":
                r = "quadrature"
            else:
                error("Invalid form representation family %s." % (form_r_family,))
        info("representation:    auto --> %s" % r)
    else:
        info("representation:    %s" % r)

    if p is None:
        p = default_precision

    # Hack to override representation with environment variable
    forced_r = os.environ.get("FFC_FORCE_REPRESENTATION")
    if forced_r:
        r = forced_r
        warning("representation:    forced by $FFC_FORCE_REPRESENTATION to '%s'" % r)
        return r, o, p

    return r, o, p


def _attach_integral_metadata(form_data, form_r_family, parameters):
    "Attach integral metadata"
    # TODO: A nicer data flow would avoid modifying the form_data at all.

    # Parameter values which make sense "per integrals" or "per integral"
    metadata_keys = (
        "representation",
        "optimize",
        # TODO: Could have finer optimize (sub)parameters here later
        "precision",
        # NOTE: We don't pass precision to quadrature, it's not
        #       worth resolving set_float_formatting hack for
        #       deprecated backends
        "quadrature_degree",
        "quadrature_rule",
    )

    # Get defaults from parameters
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

        # Determine representation, must be equal for all integrals on
        # same subdomain
        r, o, p = _determine_representation(integral_metadatas, ida, form_data,
                                            form_r_family, parameters)
        for i, integral in enumerate(ida.integrals):
            integral_metadatas[i]["representation"] = r
            integral_metadatas[i]["optimize"] = o
            integral_metadatas[i]["precision"] = p
        ida.metadata["representation"] = r
        ida.metadata["optimize"] = o
        ida.metadata["precision"] = p

        # Determine automated updates to metadata values
        for i, integral in enumerate(ida.integrals):
            qr = _autoselect_quadrature_rule(integral_metadatas[i], integral,
                                             form_data)
            qd = _autoselect_quadrature_degree(integral_metadatas[i], integral,
                                               form_data)
            integral_metadatas[i]["quadrature_rule"] = qr
            integral_metadatas[i]["quadrature_degree"] = qd

        # Extract common metadata for integral collection
        qr = _extract_common_quadrature_rule(integral_metadatas)
        qd = _extract_common_quadrature_degree(integral_metadatas)
        ida.metadata["quadrature_rule"] = qr
        ida.metadata["quadrature_degree"] = qd

        # Add common num_cells (I believe this must be equal but I'm
        # not that into this work)
        num_cells = set(md.get("num_cells") for md in integral_metadatas)
        if len(num_cells) != 1:
            error("Found integrals with different num_cells metadata on same subdomain: %s" % (str(list(num_cells)),))
        num_cells, = num_cells
        ida.metadata["num_cells"] = num_cells

        # Reconstruct integrals to avoid modifying the input integral,
        # which would affect the signature computation if the integral
        # was used again in the user program.  Modifying attributes of
        # form_data.integral_data is less problematic since it's
        # lifetime is internal to the form compiler pipeline.
        for i, integral in enumerate(ida.integrals):
            ida.integrals[i] = integral.reconstruct(metadata=integral_metadatas[i])

        # Collect all quad schemes
        quad_schemes.extend([md["quadrature_rule"] for md in integral_metadatas])

    # Validate consistency of schemes for QuadratureElements
    # TODO: Can loosen up this a bit, only needs to be consistent
    # with the integrals that the elements are used in
    _validate_quadrature_schemes_of_elements(quad_schemes,
                                             form_data.unique_sub_elements)


def _validate_quadrature_schemes_of_elements(quad_schemes, elements):
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


def _find_compatible_representations(integrals, elements):
    """Automatically select a suitable representation for integral.  Note
    that the selection is made for each integral, not for each
    term. This means that terms which are grouped by UFL into the same
    integral (if their measures are equal) will necessarily get the
    same representation.

    """
    # All representations
    compatible = set(("uflacs", "quadrature", "tsfc"))

    # Check for non-affine meshes
    if _has_higher_order_geometry(integrals):
        compatible &= set(("uflacs", "tsfc"))

    # Custom integrals
    if _has_custom_integrals(integrals):
        compatible &= set(("quadrature",))

    # Use quadrature for vertex integrals
    if any(integral.integral_type() == "vertex" for integral in integrals):
        # TODO: Test with uflacs, I think this works fine now:
        compatible &= set(("quadrature", "uflacs", "tsfc"))

    # Get ALL sub elements, needed to check for restrictions of
    # EnrichedElements.
    sub_elements = []
    for e in elements:
        sub_elements += _get_sub_elements(e)

    # Use quadrature representation if we have a quadrature element
    if any(e.family() == "Quadrature" for e in sub_elements):
        # TODO: Test with uflacs, might need a little adjustment:
        compatible &= set(("quadrature", "uflacs", "tsfc"))

    return compatible
