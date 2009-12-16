"""This module contains all the code generators that used to be in the modules:
codegenerator, finiteelement, dofmap, form and integrals. The helper functions
 (__foo_bar()) for these modules are located in codegenerators_utils"""

__author__ = "Anders Logg (logg@simula.no) and Kristian B. Oelgaard (k.b.oelgard@gmail.com)"
__date__ = "2009-12-09"
__copyright__ = "Copyright (C) 2009 Anders Logg and Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-16

# Python modules.
import numpy

# FFC modules.
from log import debug
from log import error
from utils import product
from quadratureelement import QuadratureElement
from evaluatebasis import evaluate_basis
from evaluatebasisderivatives import evaluate_basis_derivatives
from createelement import create_element

# Utility functions.
from codegenerators_utils import inner_product
from codegenerators_utils import indent
from codegenerators_utils import IndentControl

# For finite_element.
from codegenerators_utils import __extract_sub_elements
from codegenerators_utils import __generate_interpolate_vertex_values
from codegenerators_utils import __generate_evaluate_dof

# For dof_map.
from codegenerators_utils import __extract_sub_dof_maps
from codegenerators_utils import __generate_needs_mesh_entities
from codegenerators_utils import __generate_global_dimension
from codegenerators_utils import __generate_tabulate_dofs
from codegenerators_utils import __generate_tabulate_facet_dofs
from codegenerators_utils import __generate_tabulate_coordinates

#------------------------------------------------------------------------------
# From codegenerator.py
def generate_common_code(form_data, format):
    "Generate common form code according to given format."

    code = {}

    # Generate code for finite elements
    code["finite_elements"] = generate_finite_elements(form_data, format)

    # Generate code for dof maps
    code["dof_maps"] = ""#generate_dof_maps(form_data, format)

    # Generate code for form
    #debug("Generating code for form...")
    #code["form"] = generate_form(form_data, format)
    #debug("done")

    return code
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# From finiteelement.py
def generate_finite_elements(form_data, format):
    "Generate code for finite elements, including recursively nested sub elements."

    debug("Generating code for finite elements...")
    code = []

    # Iterate over form elements
    for (i, element) in enumerate(form_data.ffc_elements):

        # Extract sub elements
        sub_elements = __extract_sub_elements(element, (i,))

        # Generate code for each element
        for (label, sub_element) in sub_elements:
            code += [(label, _generate_finite_element(sub_element, format))]

    debug("done")
    return code

def _generate_finite_element(element, format):
    """Generate dictionary of code for the given finite element
    according to the given format."""

    code = {}

    # Generate code for signature
    code["signature"] = repr(element)

    # Generate code for cell_shape
    code["cell_shape"] = format["cell shape"](element.cell_domain())

    # Generate code for space_dimension
    code["space_dimension"] = "%d" % element.space_dimension()

    # Generate code for value_rank
    # FIXME: This is just a temporary hack to 'support' tensor elements
    code["value_rank"] = "%d" % element.value_dimension(0)
    #code["value_rank"] = "%d" % element._rank

    # Generate code for value_dimension
    code["value_dimension"] = ["%d" % element.value_dimension(i) for i in range(max(element.value_rank(), 1))]

    # Disable code generation for unsupported functions of QuadratureElement,
    # (or MixedElements including QuadratureElements)
    if not True in [isinstance(e, QuadratureElement) for e in element.extract_elements()]:
        # Generate code for evaluate_basis
        code["evaluate_basis"] = []#evaluate_basis(element, format)

        # Generate code for evaluate_basis_derivatives
        code["evaluate_basis_derivatives"] = []#evaluate_basis_derivatives(element, format)

        # Generate vectorised version of evaluate functions
        code["evaluate_basis_all"] =\
          format["exception"]("The vectorised version of evaluate_basis() is not yet implemented.")
        code["evaluate_basis_derivatives_all"] =\
          format["exception"]("The vectorised version of evaluate_basis_derivatives() is not yet implemented.")

        # Generate code for interpolate_vertex_values
        code["interpolate_vertex_values"] = __generate_interpolate_vertex_values(element, format)
    else:
        code["evaluate_basis"] =\
          format["exception"]("evaluate_basis() is not supported for QuadratureElement")
        code["evaluate_basis_derivatives"] =\
          format["exception"]("evaluate_basis_derivatives() is not supported for QuadratureElement")

        # Generate vectorised version of evaluate functions
        code["evaluate_basis_all"] =\
          format["exception"]("evaluate_basis_all() is not supported for QuadratureElement.")
        code["evaluate_basis_derivatives_all"] =\
          format["exception"]("evaluate_basis_derivatives_all() is not supported for QuadratureElement.")

        code["interpolate_vertex_values"] =\
          format["exception"]("interpolate_vertex_values() is not supported for QuadratureElement")

    # Generate code for evaluate_dof
    code["evaluate_dof"] = []#__generate_evaluate_dof(element, format)

    # Generate code for num_sub_elements
    code["num_sub_elements"] = "%d" % element.num_sub_elements()

    return code
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# From dofmap.py
def generate_dof_maps(form_data, format):
    "Generate code for dof maps, including recursively nested dof maps."

    debug("Generating code for finite dof maps...")
    code = []

    # Iterate over form dof maps
    for (i, dof_map) in enumerate(form_data.ffc_dof_maps):

        # Extract sub dof maps
        sub_dof_maps = __extract_sub_dof_maps(dof_map, (i,))

        # Generate code for each dof map
        for (label, sub_dof_map) in sub_dof_maps:
            code += [(label, _generate_dof_map(sub_dof_map, format))]

    debug("done")
    return code

def _generate_dof_map(dof_map, format):
    """Generate dictionary of code for the given dof map according to
    the given format"""

    code = {}

    # Generate code for signature
    code["signature"] = dof_map.signature()

    # Generate code for needs_mesh_entities
    code["needs_mesh_entities"] = __generate_needs_mesh_entities(dof_map, format)

    # Generate code for global_dimension
    code["global_dimension"] = __generate_global_dimension(dof_map, format)

    # Generate code for local_dimension
    code["local_dimension"] = "%d" % dof_map.local_dimension()

    # Generate code for geometric_dimension
    code["geometric_dimension"] = "%d" % dof_map.geometric_dimension()

    # Generate code for num_facet_dofs
    code["num_facet_dofs"] = "%d" % dof_map.num_facet_dofs()

    # Generate code for tabulate_dofs
    code["tabulate_dofs"] = __generate_tabulate_dofs(dof_map, format)

    # Generate code for tabulate_facet_dofs
    code["tabulate_facet_dofs"] = __generate_tabulate_facet_dofs(dof_map, format)

    # Generate code for tabulate_coordinates
    code["tabulate_coordinates"] = __generate_tabulate_coordinates(dof_map, format)

    # Generate code for num_sub_dof_maps
    code["num_sub_dof_maps"] = "%d" % dof_map.num_sub_dof_maps()

    return code
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# From form.py
def generate_form(form_data, format):
    """Generate dictionary of code for the given form data map
    according to the given format"""

    code = {}

    # Generate code for signature
    code["signature"] = form_data.signature

    # Generate code for rank
    code["rank"] = "%d" % form_data.rank

    # Generate code for num_coefficients
    code["num_coefficients"] = "%d" % form_data.num_coefficients

    # Generate code for num_cell_integrals
    code["num_cell_integrals"] = "%d" % form_data.num_cell_domains

    # Generate code for num_exterior_facet_integrals
    code["num_exterior_facet_integrals"] = "%d" % form_data.num_exterior_facet_domains

    # Generate code for num_interior_facet_integrals
    code["num_interior_facet_integrals"] = "%d" % form_data.num_interior_facet_domains

    return code
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# From integrals.py
def generate_combined_code(codes, form_data, prefix, format):
    "Generate combined codes for tabulation of integrals."

    combined_code = {}

    # Generate combined code for cell integrals
    args = (codes, form_data, prefix, format, "cell_integral", form_data.num_cell_domains)
    combined_code["cell_integrals"] = _generate_combined_code_common(*args)

    # Generate combined code for exterior facet integrals
    args = (codes, form_data, prefix, format, "exterior_facet_integral", form_data.num_exterior_facet_domains)
    combined_code["exterior_facet_integrals"] = _generate_combined_code_common(*args)

    # Generate combined code for interior facet integrals
    args = (codes, form_data, prefix, format, "interior_facet_integral", form_data.num_interior_facet_domains)
    combined_code["interior_facet_integrals"] = _generate_combined_code_common(*args)

    return combined_code

def _generate_combined_code_common(codes, form_data, prefix, format, integral_type, num_integrals):
    "Generate combined codes for tabulation of integrals."

    combined_code = []

    # Iterate over sub domains
    for i in range(num_integrals):

        # Add code for all representations
        contributions = []
        for (j, code) in enumerate(codes):
            key = (integral_type, i)
            if key in code:
                postfix = "%d_%s" % (i, code["representation"])
                combined_code.append((postfix, code[key]))
                contributions.append(postfix)

        # Add common code to sum up all contributions
        code = _generate_total_integral(integral_type, contributions, form_data, prefix, format)
        postfix = "%d" % i
        combined_code.append((postfix, code))

    return combined_code

def _generate_total_integral(integral_type, contributions, form_data, prefix, format):
    "Generate code for total tensor, summing contributions."

    code = {}

    # Add members
    code["members"] = ["\nprivate:\n"]
    for postfix in contributions:
        code["members"].append("  <form prefix>_%s_%s integral_%s;" % (integral_type, postfix, postfix))
    code["members"].append("")

    # FIXME: Possible optimization here not to reset entries if not needed

    # Reset all entries
    code["tabulate_tensor"] = []
    if integral_type == "interior_facet_integral":
        dims = [create_element(v.element()).space_dimension()*2 for v in form_data.arguments]
    else:
        dims = [create_element(v.element()).space_dimension() for v in form_data.arguments]
    num_entries = product(dims)
    code["tabulate_tensor"] += _generate_reset_tensor(num_entries, format)

    # Sum contributions
    code["tabulate_tensor"].append("")
    code["tabulate_tensor"].append(format["comment"]("Add all contributions to element tensor"))
    for postfix in contributions:
        # FIXME: This is UFC specific
        if integral_type == "cell_integral":
            code["tabulate_tensor"] += ["integral_%s.tabulate_tensor(A, w, c);" % postfix]
        elif integral_type == "exterior_facet_integral":
            code["tabulate_tensor"] += ["integral_%s.tabulate_tensor(A, w, c, facet);" % postfix]
        else:
            code["tabulate_tensor"] += ["integral_%s.tabulate_tensor(A, w, c0, c1, facet0, facet1);" % postfix]

    return code

def _generate_reset_tensor(num_entries, format):
    "Generate code for resetting the entries of the local element tensor."

    # Generate code as a list of declarations
    code = []

    # Comment
    code.append(format["comment"]("Reset values of the element tensor block"))

    # If number of entries is 1, just the entry equal to zero, otherwise use loop
    if num_entries == 1:
        format_element_tensor = format["element tensor"]
        format_floating_point = format["floating point"]
        code += [(format["element tensor"](0), format["floating point"](0.0))]
    else:
        # Create loop
        var = format["first free index"]
        code += [format["loop"](var, 0, num_entries)]
        name = format["element tensor quad"] + format["array access"](var)
        value = format["floating point"](0.0)
        # Use indent control to format loop code
        Indent = IndentControl()
        Indent.increase()
        code += [(Indent.indent(name), value)]
        Indent.decrease()

    return code
#------------------------------------------------------------------------------

