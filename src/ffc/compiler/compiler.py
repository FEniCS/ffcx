"""This is the compiler, taking a multi-linear form expressed as a Sum
and building the data structures (geometry and reference tensors) for
the evaluation of the multi-linear form."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-17 -- 2007-01-12"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2060
# Modified by Kristian Oelgaard 2006

# Python modules
import sys

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *
from ffc.common.exceptions import *

sys.path.append("../../")

# FFC form language modules
from ffc.formlanguage.operators import *

# FFC compiler modules
from form import *
from index import *
from tokens import *
from algebra import *
from integral import *
from signature import *
from elementsearch import *
from finiteelement import *
from dofmap import *
from mixedelement import *
from elementtensor import *
from exteriorfacettensor import *
from interiorfacettensor import *
from projection import *

def compile(sums, name = "Form", language = FFC_LANGUAGE, options = FFC_OPTIONS):
    """Compile variational form(s). This function takes as argument a
    Sum or a list of Sums representing the multilinear form(s). The
    return value is a Form or a list of Forms. Calling this function
    is equivalent to first calling build() followed by write()."""

    # Add default values for any missing options
    for key in FFC_OPTIONS:
        if not key in options:
            options[key] = FFC_OPTIONS[key]

    # Build data structures
    forms = build(sums, name, language, options)

    # Generate code
    if not forms == None:
        write(forms, options)

    return forms

def build(sums, name = "Form", language = FFC_LANGUAGE, options = FFC_OPTIONS):
    "Build data structures for evaluation of the variational form(s)."

    # Add default values for any missing options
    for key in FFC_OPTIONS:
        if not key in options:
            options[key] = FFC_OPTIONS[key]

    # Create a Form from the given sum(s)
    if isinstance(sums, list):
        forms = [Form(Sum(sum), name) for sum in sums if not sum == None]
    else:
        forms = [Form(Sum(sums), name)]

    # Check that the list is not empty
    if not len(forms) > 0:
        debug("No forms specified, nothing to do.")
        return None

    # Initialize format
    format = __choose_format(language)
    format.init(options)

    # Generate the element tensor for all given forms
    for form in forms:
        __build_form(form, format, options)

    # Return form(s)
    if len(forms) > 1:
        return forms
    else:
        return forms[0]

def write(forms, options = FFC_OPTIONS):
    "Generate code from previously built data structures."

    # Add default values for any missing options
    for key in FFC_OPTIONS:
        if not key in options:
            options[key] = FFC_OPTIONS[key]

    # Make sure we have a list of forms
    if isinstance(forms, list):
        forms = [form for form in forms if not form == None]
    elif not forms == None:
        forms = [forms]
    else:
        forms = []

    # Check that the list is not empty
    if not len(forms) > 0:
        debug("No forms specified, nothing to do.", 1)
        return None

    # Get output format (all forms have the same format, so pick the first)
    format = forms[0].format

    # Generate output for all forms
    format.write(forms, options)

    return

def writeFiniteElement(element, name = "MyElement", language = FFC_LANGUAGE, options = FFC_OPTIONS):
    "Generate code for the given finite element."

    # Add default values for any missing options
    for key in FFC_OPTIONS:
        if not key in options:
            options[key] = FFC_OPTIONS[key]

    # Initialize format
    format = __choose_format(language)
    format.init(options)

    # Write a nice message
    debug("Compiling finite element: " + str(element))

    # Generate code
    format.writeFiniteElement(element, name, options)

def __build_form(form, format, options):
    "Build data structures for evaluation of the variational form."
    
    debug("\nCompiling form: " + str(form), 0)
    debug("Number of terms in form: %d" % len(form.sum.products), 1)
        
    # Count the number of functions
    form.nfunctions = max_index(form.sum, Index.FUNCTION) + 1
    debug("Number of functions (coefficients): " + str(form.nfunctions), 1)

    # Count the number of projections
    form.nprojections = max_index(form.sum, Index.PROJECTION) + 1
    debug("Number of projections (coefficients): " + str(form.nprojections), 1)

    # Count the number of constants
    form.nconstants = max_index(form.sum, Index.CONSTANT) + 1
    debug("Number of constants: " + str(form.nconstants), 1)

    # Find the test and trial finite elements
    form.test = find_test(form.sum)
    form.trial = find_trial(form.sum)

    # Find the original elements for all functions
    form.elements = find_elements(form.sum, form.nfunctions)

    # Find the projections for all functions
    form.projections = find_projections(form.sum, form.nprojections)

    # Create empty sets of used coefficient declarations
    cK_used  = set()
    cSe_used = set()
    cSi_used = set()

    # Create empty sets of used geometry tensor declarations
    gK_used  = set()
    gSe_used = set()
    gSi_used = set()

    # Compute element tensor for cell
    debug("Compiling tensor representation over cells")
    form.AK = ElementTensor(form.sum, format, cK_used, gK_used, options)
    form.num_ops = form.AK.num_ops
        
    # FIXME: Number of operations not counted for facet terms

    # Compute element tensors for exterior facets
    debug("Compiling tensor representation over exterior facets")
    num_facets = form.sum.products[0].basisfunctions[0].element.num_facets()
    form.ASe = [None for i in range(num_facets)]
    for i in range(num_facets):
        form.ASe[i] = ExteriorFacetTensor(form.sum, format, cSe_used, gSe_used, options, i)

    # Compute element tensors for combinations of interior facets
    debug("Compiling tensor representation over exterior facets")
    num_facets = form.sum.products[0].basisfunctions[0].element.num_facets()
    num_alignments = form.sum.products[0].basisfunctions[0].element.num_alignments()
    form.ASi = [[[None for k in range(num_alignments)] for j in range(num_facets)] for i in range(num_facets)]
    for i in range(num_facets):
        for j in range (num_facets):
            for k in range(num_alignments):
                form.ASi[i][j][k] = InteriorFacetTensor(form.sum, format, cSi_used, gSi_used, options, i, j, k)

    # Report number of operations
    debug("Number of operations (multiplications) in computation of element tensor: " + str(form.num_ops), 1)

    # Compute coefficient declarations, common to all terms
    form.cK  = __compute_coefficients(form.projections, format, cK_used,  1)
    form.cSe = __compute_coefficients(form.projections, format, cSe_used, 1)
    form.cSi = __compute_coefficients(form.projections, format, cSi_used, 2)
      
    # Check primary ranks
    __check_primary_ranks(form, num_facets, num_alignments)
    
    # Save format
    form.format = format

    # Copy variables to new names
    # FIXME: Move to new names and put this code into the Form class
    form.signature = str(form.sum) # FIXME: Use signature.py
    form.rank = form.rank
    form.num_coefficients = form.nfunctions
    form.finite_elements = []
    if not form.test == None:
        form.finite_elements += [form.test]
    if not form.trial == None:
        form.finite_elements += [form.trial]
    form.finite_elements += form.elements
    form.dof_maps = []
    from ffc.format import ufcformat
    if format == ufcformat:
        for element in form.finite_elements:
            form.dof_maps += [DofMap(element, format)]

def __check_primary_ranks(form, num_facets, num_alignments):
    "Check that all primary ranks are equal."

    form.rank = None
    form.dims = None
    form.indices = None

    # Check ranks for cell terms
    for term in form.AK.terms:
        if form.rank == None:
            form.rank = term.A0.i.rank
            form.dims = term.A0.i.dims
            form.indices = term.A0.i.indices
        elif not form.rank == term.A0.i.rank:
            raise FormError(form.sum, "Form must be linear in each of its arguments.")

    # Check ranks for exterior facet terms
    for i in range(num_facets):
        for term in form.ASe[i].terms:
            if form.rank == None:
                form.rank = term.A0.i.rank
                form.dims = term.A0.i.dims
                form.indices = term.A0.i.indices
            elif not form.rank == term.A0.i.rank:
                raise FormError(form.sum, "Form must be linear in each of its arguments.")

    # Check ranks for interior facet terms
    for i in range(num_facets):
        for j in range(num_facets):
            for k in range(num_alignments):
                for term in form.ASi[i][j][k].terms:
                    if form.rank == None:
                        form.rank = term.A0.i.rank
                        form.dims = term.A0.i.dims
                        form.indices = term.A0.i.indices
                    elif not form.rank == term.A0.i.rank:
                        raise FormError(form.sum, "Form must be linear in each of its arguments.")

def __compute_coefficients(projections, format, c_used, num_spaces):
    "Precompute declarations of coefficients according to given format."

    declarations = []

    # Iterate over projections
    for (n0, n1, e0, e1, P) in projections:

        if P == None:
            n = num_spaces*e0.space_dimension()
            # No projection, just copy the values
            for k in range(n):
                name = format.format["coefficient"](n1, k)
                value = format.format["coefficient table"](n0, k)
                declaration = Declaration(name, value)
                # Mark entries that are used
                if name in c_used:
                    declaration.used = True
                declarations += [declaration]
        else:
            for i in range(num_spaces):
                # Compute projection
                (m, n) = numpy.shape(P)
                for k in range(m):
                    terms = []
                    for l in range(n):
                        if abs(P[k][l] < FFC_EPSILON):
                            continue
                        cl = format.format["coefficient table"](n0, l + i*n)
                        if abs(P[k][l] - 1.0) < FFC_EPSILON:
                            terms += [cl]
                        else:
                            Pkl = format.format["floating point"](P[k][l])
                            terms += [format.format["multiplication"]([Pkl, cl])]
                    name = format.format["coefficient"](n1, k + i*m)
                    value = format.format["sum"](terms)
                    declaration = Declaration(name, value)
                    # Mark entries that are used
                    if name in c_used:
                        declaration.used = True
                    # Add to list of declarations
                    declarations += [declaration]

    return declarations

def __choose_format(language):
    "Choose format from specified language."

    # Default format
    if not language:
        language = "dolfin"

    # Choose format
    if language.lower() == "ufc":
        from ffc.format import ufcformat
        format = ufcformat
    elif language.lower() == "dolfin":
        from ffc.format import dolfin
        format = dolfin
    elif language.lower() == "latex":
        from ffc.format import latex
        format = latex
    elif language.lower() == "raw":
        from ffc.format import raw
        format = raw
    elif language.lower() == "ase":
        from ffc.format import ase
        format = ase
    elif language.lower() == "xml":
        from ffc.format import xml
        format = xml
    else:
        raise "RuntimeError", "Unknown language " + str(language)
    return format
