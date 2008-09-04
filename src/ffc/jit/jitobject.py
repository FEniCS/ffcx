__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2008-09-04 -- 2008-09-04"
__copyright__ = "Copyright (C) 2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from hashlib import sha1
from ffc.compiler.analysis import simplify, analyze

# Used for caching
_jit_objects = {}

class JITObject:
    """This class is a wrapper for compiled objects in the context of
    specific compiler options. JITObjects with the same signature are
    considered to be identical and may be mapped to the same code.
    This allows code and modules to be reused from cache."""

    def __init__(self, form, options):
        "Create JITObject for given form and options"
        self.form = form
        self.options = options
        self._form_data = None
        self._signature = None

    def signature(self):
        "Return signature"

        # Check if we have computed the signature before
        if not self._signature is None:
            return self._signature

        # Compute signature
        self.form_data = analyze.analyze(self.form, simplify_form=False)
        form_signature = str(self.form)
        element_signature = ";".join([element.signature() for element in self.form_data.elements])
        options_signature = str(self.options)
        s = ";".join([form_signature, element_signature, options_signature])
        self._signature = "form_" + sha1(s).hexdigest()
        
        return self._signature

def wrap(form, options):
    "Wrap form into a JITObject"
    key = (form, str(options))
    if key in _jit_objects:
        jit_object = _jit_objects[key]
    else:
        jit_object = JITObject(form, options)
        _jit_objects[key] = jit_object
    return jit_object
