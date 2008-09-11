__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2008-09-04 -- 2008-09-04"
__copyright__ = "Copyright (C) 2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from hashlib import sha1

# FFC compiler modules
from ffc.compiler.language import algebra
from ffc.compiler.analysis import simplify, analyze

class JITObject:
    """This class is a wrapper for a compiled object in the context of
    specific compiler options. A JITObject is identified either by its
    hash value or by its signature. The hash value is valid only in a
    single instance of an application (at runtime). The signature is
    persistent and may be used for caching modules on disk."""

    def __init__(self, form, options):
        "Create JITObject for given form and options"
        self.form = form
        self.options = options
        self._form_data = None
        self._hash = None
        self._signature = None

    def __hash__(self):
        "Return unique integer for form + options"

        # Check if we have computed the hash before
        if not self._hash is None:
            return self._hash

        # Compute hash (use a hack to convert hexdigest to int)
        string = str(id(self.form)) + str(self.options)
        hexdigest = sha1(string).hexdigest()
        number = int("".join([c for c in hexdigest if c.isdigit()]))
        
        return number

    def __eq__(self, other):
        "Check for equality"
        return hash(self) == hash(other)

    def signature(self):
        "Return unique string for form expression + options"

        # Check if we have computed the signature before
        if not self._signature is None:
            return self._signature

        # Compute signature
        self.form_data = analyze.analyze(algebra.Form(self.form), simplify_form=False)
        form_signature = str(self.form)
        element_signature = ";".join([element.signature() for element in self.form_data.elements])
        options_signature = str(self.options)
        string = ";".join([form_signature, element_signature, options_signature])
        self._signature = "form_" + sha1(string).hexdigest()
        
        return self._signature
