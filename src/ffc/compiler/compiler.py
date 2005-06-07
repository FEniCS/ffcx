"""This is the compiler, taking a multi-linear form expressed as a Sum
and building the data structures (geometry and reference tensors) for
the evaluation of the multi-linear form."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-17 -- 2005-05-20"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys
from Numeric import *

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC format modules
sys.path.append("../../")
from ffc.format import dolfin
from ffc.format import latex
from ffc.format import raw

# FFC compiler modules
from form import *
from index import *
from algebra import *
from integral import *
from reassign import *
from elementsearch import *
from finiteelement import *
from elementtensor import *

def compile(sums, name = "Form", language = "C++", license = FFC_LICENSE):
    """Generate code for evaluation of the variational form.
    This function takes as argument a Sum or a list of Sums
    representing the multilinear form(s). The return value is a Form
    or a list of Forms."""

    # Create a Form from the given sum(s)
    if isinstance(sums, list):
        forms = [Form(Sum(sum), name) for sum in sums if sum]
    else:
        forms = [Form(Sum(sums), name)]

    # Choose language
    if not language:
        format = dolfin
    elif language == "C++" or language == "c++":
        format = dolfin
    elif language == "LaTeX" or language == "latex":
        format = latex
    elif language == "raw":
        format = raw
    elif language == "ase":
        from ffc.format import ase
        format = ase
    else:
        raise "RuntimeError", "Unknown language " + str(language)

    # Generate the element tensor for all given forms
    for form in forms:

        debug("\nCompiling form: " + str(form), 0)
        debug("Number of terms in form: %d" % len(form.sum.products), 1)
        
        # Count the number of functions
        form.nfunctions = max_index(form.sum, "function") + 1
        debug("Number of functions (coefficients): " + str(form.nfunctions), 1)

        # Count the number of constants
        form.nconstants = max_index(form.sum, "constant") + 1
        debug("Number of constants: " + str(form.nconstants), 1)

        # Find the test, trial and function finite elements
        form.test = find_test(form.sum)
        form.trial = find_trial(form.sum)
        #(form.elements, form.felement) = find_elements(form.sum, form.nfunctions)
        form.elements = find_elements(form.sum, form.nfunctions)

        # Create element tensors
        form.AKi = ElementTensor(form.sum, "interior", format)
        form.AKb = ElementTensor(form.sum, "boundary", format)

        # Check primary ranks
        __check_primary_ranks(form)

    # Generate output
    format.compile(forms, license)

    # Return form
    if len(forms) > 1:
        forms
    else:
        return forms[0]

def __check_primary_ranks(form):
    "Check that all primary ranks are equal."
    terms = form.AKi.terms + form.AKb.terms
    ranks = [term.A0.i.rank for term in terms]
    if not ranks[1:] == ranks[:-1]:
        "Form must be linear in each of its arguments."
    form.rank = ranks[0]
    form.dims = terms[0].A0.i.dims
    form.indices = terms[0].A0.i.indices
    return

if __name__ == "__main__":

    print "Testing form compiler"
    print "---------------------"

    element = FiniteElement("Lagrange", "triangle", 1)
    
    u = BasisFunction(element)
    v = BasisFunction(element)
    i = Index()
    dx = Integral("interior")
    ds = Integral("boundary")
    
    a = u.dx(i)*v.dx(i)*dx + u*v*ds
    compile(a, "form", "C++")
    compile(a, "form", "LaTeX")
    compile(a, "form", "raw")
