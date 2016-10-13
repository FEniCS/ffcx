# -*- coding: utf-8 -*-
# Copyright (C) 2015-2016 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

#import inspect
import re
#from string import Formatter

#from ufl import product
from ffc.log import error, warning
#from ffc.backends.ufc import *

import ffc.backends.ufc

from ffc.uflacs.language.format_lines import format_indented_lines
from ffc.uflacs.backends.ufc.templates import *

#__all__ = (["ufc_form", "ufc_dofmap", "ufc_finite_element", "ufc_integral"]
#           + ["ufc_%s_integral" % integral_type for integral_type in integral_types])


# These are all the integral types directly supported in ufc
from ffc.representation import ufc_integral_types


# These are the method names in ufc::form that are specialized for each integral type
integral_name_templates = (
    "max_%s_subdomain_id",
    "has_%s_integrals",
    "create_%s_integral",
    "create_default_%s_integral",
    )


class ufc_generator(object):
    """Common functionality for code generators producing ufc classes.

    The generate function is the driver for generating code for a class.
    It automatically extracts template keywords and inserts the results
    from calls to self.<keyword>(language, ir), or the value of ir[keyword]
    if there is no self.<keyword>.
    """
    def __init__(self, basename):
        ufc_templates = ffc.backends.ufc.templates
        self._header_template = ufc_templates[basename + "_header"]
        self._implementation_template = ufc_templates[basename + "_implementation"]
        self._combined_template = ufc_templates[basename + "_combined"]

        r = re.compile(r"%\(([a-zA-Z0-9_]*)\)")
        self._header_keywords = set(r.findall(self._header_template))
        self._implementation_keywords = set(r.findall(self._implementation_template))
        self._combined_keywords = set(r.findall(self._combined_template))

        self._keywords = sorted(self._header_keywords | self._implementation_keywords)

        # Do some ufc interface template checking, to catch bugs early when we change the ufc interface templates
        if set(self._keywords) != set(self._combined_keywords):
            a = set(self._header_keywords) - set(self._combined_keywords)
            b = set(self._implementation_keywords) - set(self._combined_keywords)
            c = set(self._combined_keywords) - set(self._keywords)
            msg = "Templates do not have matching sets of keywords:"
            if a:
                msg += "\n  Header template keywords '%s' are not in the combined template." % (sorted(a),)
            if b:
                msg += "\n  Implementation template keywords '%s' are not in the combined template." % (sorted(b),)
            if c:
                msg += "\n  Combined template keywords '%s' are not in the header or implementation templates." % (sorted(c),)
            error(msg)

    def generate_snippets(self, L, ir):
        "Generate code snippets for each keyword found in templates."
        snippets = {}
        for kw in self._keywords:
            # Check that attribute self.<keyword> is available
            if not hasattr(self, kw):
                error("Missing handler for keyword '%s' in class %s." % (kw, self.__class__.__name__))

            # Call self.<keyword>(L, ir) to get value
            method = getattr(self, kw)
            value = method(L, ir)

            # Indent body and format to str
            if isinstance(value, L.CStatement):
                value = L.Indented(value.cs_format())
                value = format_indented_lines(value)

            # Store formatted code in snippets dict
            snippets[kw] = value

        # Error checking (can detect some bugs early when changing the ufc interface)
        # Get all attributes of subclass class (skip "_foo")
        attrs = set(name for name in dir(self) if not (name.startswith("_") or name.startswith("generate")))
        # Get all attributes of this base class (skip "_foo" and "generate*")
        base_attrs = set(name for name in dir(ufc_generator) if not (name.startswith("_") or name.startswith("generate")))
        # The template keywords should not contain any names not among the class attributes
        missing = set(self._keywords) - attrs
        if missing:
            warning("*** Missing generator functions:\n%s" % ('\n'.join(map(str, sorted(missing))),))
        # The class attributes should not contain any names not among the template keywords
        # (this is strict, a useful check when changing ufc, but can be dropped)
        unused = attrs - set(self._keywords) - base_attrs
        if unused:
            warning("*** Unused generator functions:\n%s" % ('\n'.join(map(str, sorted(unused))),))

        # Return snippets, a dict of code strings
        return snippets

    def generate(self, L, ir, snippets=None):
        "Return composition of templates with generated snippets."
        if snippets is None:
            snippets = self.generate_snippets(L, ir)
        h = self._header_template % snippets
        cpp = self._implementation_template % snippets
        return h, cpp

    def classname(self, L, ir):
        "Return classname."
        return ir["classname"]

    def members(self, L, ir):
        "Return empty string. Override in classes that need members."
        if ir.get("members"):
            error("Missing generator function.")
        return ""

    def constructor(self, L, ir):
        "Return empty string. Override in classes that need constructor."
        if ir.get("constructor"):
            error("Missing generator function.")
        return ""

    def constructor_arguments(self, L, ir):
        "Return empty string. Override in classes that need constructor."
        if ir.get("constructor_arguments"):
            error("Missing generator function.")
        return ""

    def initializer_list(self, L, ir):
        "Return empty string. Override in classes that need constructor."
        if ir.get("initializer_list"):
            error("Missing generator function.")
        return ""

    def destructor(self, L, ir):
        "Return empty string. Override in classes that need destructor."
        if ir.get("destructor"):
            error("Missing generator function.")
        return ""

    def signature(self, L, ir):
        "Default implementation of returning signature string fetched from ir."
        sig = ir["signature"]
        return L.Return(L.LiteralString(sig))

    def create(self, L, ir):
        "Default implementation of creating a new object of the same type."
        classname = ir["classname"]
        return L.Return(L.New(classname))
