"""This is the compiler, taking a multi-linear form expressed as a Form
and building the data structures (geometry and reference tensors) for
the evaluation of the multi-linear form."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-17 -- 2007-01-23"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2060
# Modified by Kristian Oelgaard 2006

# Python modules
import sys
from sets import Set

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *
from ffc.common.exceptions import *

sys.path.append("../../")

# FFC form language modules
from ffc.formlang.index import *
from ffc.formlang.tokens import *
from ffc.formlang.algebra import *
from ffc.formlang.integral import *
from ffc.formlang.signature import *
from ffc.formlang.operators import *

# FFC codegen modules
from ffc.codegen.codegenerator import *

# FFC format modules
from ffc.format import dolfin
from ffc.format import ufcformat

# FFC fem modules
from ffc.fem.mixedelement import *
from ffc.fem.dofmap import *

# FFC compiler modules
from form import *
from elementsearch import *
from elementtensor import *
from exteriorfacettensor import *
from interiorfacettensor import *
from projection import *

def compile(forms, name = "Form", language = FFC_LANGUAGE, options = FFC_OPTIONS):
    """Compile variational form(s). This function takes as argument a
    Form or a list of Forms representing the multilinear form(s). The
    return value is a FormCode or a list of FormCodes. Calling this
    function is equivalent to first calling build() followed by
    write()."""

    # Add default values for any missing options
    for key in FFC_OPTIONS:
        if not key in options:
            options[key] = FFC_OPTIONS[key]

    # Build data structures
    (formcodes, code) = build(forms, name, language, options)

    # Generate code
    if not formcodes == None:
        write(formcodes, code, options)

    return formcodes

def build(forms, name = "Form", language = FFC_LANGUAGE, options = FFC_OPTIONS):
    "Build data structures for evaluation of the variational form(s)."

    # Add default values for any missing options
    for key in FFC_OPTIONS:
        if not key in options:
            options[key] = FFC_OPTIONS[key]

    # Create a FormCode from the given Form(s)
    if isinstance(forms, list):
        formcodes = [FormCode(Form(form), name) for form in forms if not form == None]
    else:
        formcodes = [FormCode(Form(forms), name)]

    # Check that the list is not empty
    if not len(formcodes) > 0:
        debug("No forms specified, nothing to do.")
        return None

    # Initialize format
    format = __choose_format(language)
    format.init(options)

    # Generate the element tensor for all given forms
    for formcode in formcodes:
        code = __build_form(formcode, format, options)

    # Return formcode(s)
    if len(formcodes) > 1:
        return (formcodes, code)
    else:
        return (formcodes[0], code)

def write(formcodes, code, options = FFC_OPTIONS):
    "Generate code from previously built data structures."

    # Add default values for any missing options
    for key in FFC_OPTIONS:
        if not key in options:
            options[key] = FFC_OPTIONS[key]

    # Make sure we have a list of formcodes
    if isinstance(formcodes, list):
        formcodes = [formcode for formcode in formcodes if not formcode == None]
    elif not formcodes == None:
        formcodes = [formcodes]
    else:
        formcodes = []

    # Check that the list is not empty
    if not len(formcodes) > 0:
        debug("No forms specified, nothing to do.", 1)
        return None

    # Get output format (all forms have the same format, so pick the first)
    format = formcodes[0].format

    # Generate output for all formcodes
    # FIXME: Should be two arguments when this is fixed: code, options
    format.write(formcodes, code, options)

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

def __build_form(formcode, format, options):
    "Build data structures for evaluation of the variational form."
    
    debug("\nCompiling form: " + str(formcode), 0)
    debug("Number of terms in form: %d" % len(formcode.form.monomials), 1)
        
    # Count the number of functions
    formcode.nfunctions = max_index(formcode.form, Index.FUNCTION) + 1
    debug("Number of functions (coefficients): " + str(formcode.nfunctions), 1)

    # Count the number of projections
    formcode.nprojections = max_index(formcode.form, Index.PROJECTION) + 1
    debug("Number of projections (coefficients): " + str(formcode.nprojections), 1)

    # Count the number of constants
    formcode.nconstants = max_index(formcode.form, Index.CONSTANT) + 1
    debug("Number of constants: " + str(formcode.nconstants), 1)

    # Find the test and trial finite elements
    formcode.test = find_test(formcode.form)
    formcode.trial = find_trial(formcode.form)

    # Find the original elements for all functions
    formcode.elements = find_elements(formcode.form, formcode.nfunctions)

    # Find the projections for all functions
    formcode.projections = find_projections(formcode.form, formcode.nprojections)

    # Create empty sets of used coefficient declarations
    cK_used  = Set()
    cSe_used = Set()
    cSi_used = Set()

    # Create empty sets of used geometry tensor declarations
    gK_used  = Set()
    gSe_used = Set()
    gSi_used = Set()

    # Compute element tensor for cell
    debug("Compiling tensor representation over cells")
    formcode.AK = ElementTensor(formcode.form, format, cK_used, gK_used, options)
    formcode.num_ops = formcode.AK.num_ops
        
    # FIXME: Number of operations not counted for facet terms

    # Compute element tensors for exterior facets
    debug("Compiling tensor representation over exterior facets")
    num_facets = formcode.form.monomials[0].basisfunctions[0].element.num_facets()
    formcode.ASe = [None for i in range(num_facets)]
    for i in range(num_facets):
        formcode.ASe[i] = ExteriorFacetTensor(formcode.form, format, cSe_used, gSe_used, options, i)

    # Compute element tensors for combinations of interior facets
    debug("Compiling tensor representation over exterior facets")
    num_facets = formcode.form.monomials[0].basisfunctions[0].element.num_facets()
    num_alignments = formcode.form.monomials[0].basisfunctions[0].element.num_alignments()
    formcode.ASi = [[[None for k in range(num_alignments)] for j in range(num_facets)] for i in range(num_facets)]
    for i in range(num_facets):
        for j in range (num_facets):
            for k in range(num_alignments):
                formcode.ASi[i][j][k] = InteriorFacetTensor(formcode.form, format, cSi_used, gSi_used, options, i, j, k)

    # Report number of operations
    debug("Number of operations (multiplications) in computation of element tensor: " + str(formcode.num_ops), 1)

    # Compute coefficient declarations, common to all terms
    formcode.cK  = __compute_coefficients(formcode.projections, format, cK_used,  1)
    formcode.cSe = __compute_coefficients(formcode.projections, format, cSe_used, 1)
    formcode.cSi = __compute_coefficients(formcode.projections, format, cSi_used, 2)
      
    # Check primary ranks
    __check_primary_ranks(formcode, num_facets, num_alignments)
    
    # Save format
    formcode.format = format

    # Copy variables to new names
    # FIXME: Move to new names and put this code into the FormCode class
    formcode.signature = str(formcode.form) # FIXME: Use signature.py
    formcode.rank = formcode.rank
    formcode.num_coefficients = formcode.nfunctions
    formcode.num_arguments = formcode.rank + formcode.num_coefficients
    formcode.finite_elements = []
    if not formcode.test == None:
        formcode.finite_elements += [formcode.test]
    if not formcode.trial == None:
        formcode.finite_elements += [formcode.trial]
    formcode.finite_elements += formcode.elements
    formcode.dof_maps = [DofMap(element) for element in formcode.finite_elements]

    # Generate code
    code = {}
    if format == ufcformat:
        code = generate_code(formcode.finite_elements, formcode.dof_maps, format)

    return code

def __check_primary_ranks(formcode, num_facets, num_alignments):
    "Check that all primary ranks are equal."

    formcode.rank = None
    formcode.dims = None
    formcode.indices = None

    # Check ranks for cell terms
    for term in formcode.AK.terms:
        if formcode.rank == None:
            formcode.rank = term.A0.i.rank
            formcode.dims = term.A0.i.dims
            formcode.indices = term.A0.i.indices
        elif not formcode.rank == term.A0.i.rank:
            raise FormcodeError(formcode.form, "Formcode must be linear in each of its arguments.")

    # Check ranks for exterior facet terms
    for i in range(num_facets):
        for term in formcode.ASe[i].terms:
            if formcode.rank == None:
                formcode.rank = term.A0.i.rank
                formcode.dims = term.A0.i.dims
                formcode.indices = term.A0.i.indices
            elif not formcode.rank == term.A0.i.rank:
                raise FormcodeError(formcode.form, "Formcode must be linear in each of its arguments.")

    # Check ranks for interior facet terms
    for i in range(num_facets):
        for j in range(num_facets):
            for k in range(num_alignments):
                for term in formcode.ASi[i][j][k].terms:
                    if formcode.rank == None:
                        formcode.rank = term.A0.i.rank
                        formcode.dims = term.A0.i.dims
                        formcode.indices = term.A0.i.indices
                    elif not formcode.rank == term.A0.i.rank:
                        raise FormcodeError(formcode.form, "Formcode must be linear in each of its arguments.")

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
        format = ufcformat
    elif language.lower() == "dolfin":
        format = dolfin
    elif language.lower() == "latex":
        format = latex
    elif language.lower() == "raw":
        format = raw
    elif language.lower() == "ase":
        format = ase
    elif language.lower() == "xml":
        format = xml
    else:
        raise "RuntimeError", "Unknown language " + str(language)
    return format
