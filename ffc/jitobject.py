# Copyright (C) 2008-2013 Anders Logg
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
# Modified by Martin Alnaes, 2013

# Python modules.
from hashlib import sha1

# UFL modules.
import ufl
from ufl.utils.sorting import canonicalize_metadata

# FFC modules.
from ffc import __version__ as FFC_VERSION
from ffc.ufc_signature import ufc_signature
from ffc.parameters import compilation_relevant_parameters

# UFC modules.
from ffc.backends import ufc


class JITObject:
    """This class is a wrapper for a compiled object in the context of
    specific compiler parameters. A JITObject is identified either by its
    hash value or by its signature. The hash value is valid only in a
    single instance of an application (at runtime). The signature is
    persistent and may be used for caching modules on disk."""

    def __init__(self, form, parameters):
        "Create JITObject for given form and parameters"
        assert isinstance(form, ufl.Form) or isinstance(form, ufl.FiniteElementBase)

        # Store data
        self.form = form
        self.parameters = parameters
        self._hash = None
        self._signature = None

    def __hash__(self):
        "Return unique integer for form + parameters"
        # Check if we have computed the hash before
        if self._hash is None:
            # Compute hash based on signature
            self._hash = int(self.signature(), 16)
        return self._hash

    def __eq__(self, other):
        "Check for equality"
        return hash(self) == hash(other)

    def signature(self):
        "Return unique string for form + parameters"

        # Check if we have computed the signature before
        if not self._signature is None:
            return self._signature

        # Get signature from form
        if isinstance(self.form, ufl.Form):
            form_signature = self.form.signature()
        elif isinstance(self.form, ufl.FiniteElementBase):
            form_signature = repr(self.form)

        # Compute other relevant signatures
        parameters_signature = _parameters_signature(self.parameters)
        ffc_signature = str(FFC_VERSION)

        # Build common signature
        signatures = [form_signature,
                      parameters_signature,
                      ffc_signature,
                      ufc_signature()]
        string = ";".join(signatures)

        self._signature = sha1(string.encode('utf-8')).hexdigest()

        # Uncomment for debugging
        #print "form_signature       =", form_signature
        #print "parameters_signature =", parameters_signature
        #print "ffc_signature        =", ffc_signature
        #print "signature            =", self._signature

        return self._signature

def _parameters_signature(parameters):
    "Return parameters signature (some parameters must be ignored)."
    parameters = compilation_relevant_parameters(parameters)
    return str(canonicalize_metadata(parameters))
