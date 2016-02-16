# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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

"""FIXME: Translate these classes to the CNode hierarchy."""


class TemplateArgumentList(ASTNode):
    singlelineseparators = ('<', ', ', '>')
    multilineseparators = ('<\n', ',\n', '\n>')

    def __init__(self, args, multiline=True):
        self.args = args
        self.multiline = multiline

    def format(self, level):
        if self.multiline:
            container = Indented
            start, sep, end = self.multilineseparators
        else:
            container = tuple
            start, sep, end = self.singlelineseparators
            # Add space to avoid >> template issue
            last = self.args[-1]
            if isinstance(last, TemplateArgumentList) or (
                    isinstance(last, Type) and last.template_arguments):
                end = ' ' + end
        code = [sep.join(format_code(arg) for arg in self.args)]
        code = (start, container(code), end)
        return format_code(code, level)


class Type(ASTNode):
    def __init__(self, name, template_arguments=None, multiline=False):
        self.name = name
        self.template_arguments = template_arguments
        self.multiline = multiline

    def format(self, level):
        code = self.name
        if self.template_arguments:
            code = code, TemplateArgumentList(self.template_arguments, self.multiline)
        return format_code(code, level)


class TypeDef(ASTNode):
    def __init__(self, type_, typedef):
        self.type_ = type_
        self.typedef = typedef

    def format(self, level):
        code = ('typedef ', self.type_, " %s;" % self.typedef)
        return format_code(code, level)


# TODO: Add variable access type with type checking to replace explicit str instances all over the place.


class ArrayAccess(ASTOperator):
    def __init__(self, arraydecl, indices):
        if isinstance(arraydecl, ArrayDecl):
            self.arrayname = arraydecl.name
        else:
            self.arrayname = arraydecl

        if isinstance(indices, (list, tuple)):
            self.indices = indices
        else:
            self.indices = (indices,)

        # Early error checking of array dimensions
        if any(isinstance(i, int) and i < 0 for i in self.indices):
            raise ValueError("Index value < 0.")

        # Additional checks possible if we get an ArrayDecl instead of just a name
        if isinstance(arraydecl, ArrayDecl):
            if len(self.indices) != len(arraydecl.sizes):
                raise ValueError("Invalid number of indices.")
            if any((isinstance(i, int) and isinstance(d, int) and i >= d)
                   for i, d in zip(self.indices, arraydecl.sizes)):
                raise ValueError("Index value >= array dimension.")

    def format(self, level):
        brackets = tuple(("[", n, "]") for n in self.indices)
        code = (self.arrayname, brackets)
        return format_code(code, level)


class Class(ASTStatement):
    def __init__(self, name, superclass=None, public_body=None,
                 protected_body=None, private_body=None,
                 template_arguments=None, template_multiline=False):
        self.name = name
        self.superclass = superclass
        self.public_body = public_body
        self.protected_body = protected_body
        self.private_body = private_body
        self.template_arguments = template_arguments
        self.template_multiline = template_multiline

    def format(self, level):
        code = []
        if self.template_arguments:
            code += [('template', TemplateArgumentList(self.template_arguments,
                                                       self.template_multiline))]
        if self.superclass:
            code += ['class %s: public %s' % (self.name, self.superclass)]
        else:
            code += ['class %s' % self.name]
        code += ['{']
        if self.public_body:
            code += ['public:', Indented(self.public_body)]
        if self.protected_body:
            code += ['protected:', Indented(self.protected_body)]
        if self.private_body:
            code += ['private:', Indented(self.private_body)]
        code += ['};']
        return format_code(code, level)
