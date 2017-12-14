# -*- coding: utf-8 -*-
# Copyright (C) 2015-2017 Martin Sandve Alnæs
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

# TODO: Move these to uflacs.language utils?


def generate_return_new(L, classname, factory):
    if factory:
        return L.Return(L.Call("create_" + classname))
    else:
        return L.Return(L.New(classname))


def generate_return_new_switch(L, i, classnames, args=None, factory=False):
    # TODO: UFC functions of this type could be replaced with return vector<shared_ptr<T>>{objects}.

    if isinstance(i, str):
        i = L.Symbol(i)

    if factory:
        def create(classname):
            return L.Call("create_" + classname)
    else:
        def create(classname):
            return L.New(classname)

    default = L.Return(L.Null())
    if classnames:
        cases = []
        if args is None:
            args = list(range(len(classnames)))
        for j, classname in zip(args, classnames):
            if classname:
                cases.append((j, L.Return(create(classname))))
        return L.Switch(i, cases, default=default)
    else:
        return default


def generate_return_literal_switch(L, i, values, default, literal_type, typename=None):
    # TODO: UFC functions of this type could be replaced with return vector<T>{values}.

    if isinstance(i, str):
        i = L.Symbol(i)
    return_default = L.Return(literal_type(default))

    if values and typename is not None:
        # Store values in static table and return from there
        V = L.Symbol("return_values")
        decl = L.ArrayDecl("static const %s" % typename, V, len(values),
                            [literal_type(k) for k in values])
        return L.StatementList([
            decl,
            L.If(L.GE(i, len(values)),
                 return_default),
            L.Return(V[i])
            ])
    elif values:
        # Need typename to create static array, fallback to switch
        cases = [(j, L.Return(literal_type(k)))
                 for j, k in enumerate(values)]
        return L.Switch(i, cases, default=return_default)
    else:
        # No values, just return default
        return return_default


def generate_return_sizet_switch(L, i, values, default):
    return generate_return_literal_switch(L, i, values, default, L.LiteralInt, "std::size_t")


def generate_return_int_switch(L, i, values, default):
    return generate_return_literal_switch(L, i, values, default, L.LiteralInt, "int")


def generate_return_bool_switch(L, i, values, default):
    return generate_return_literal_switch(L, i, values, default, L.LiteralBool, "bool")


# TODO: Better error handling
def generate_error(L, msg, emit_warning):
    if emit_warning:
        return L.VerbatimStatement('std::cerr << "*** FFC warning: " << "%s" << std::endl;' % (msg,))
    else:
        return L.Throw('std::runtime_error', msg)


"""
// TODO: UFC functions that return constant arrays
// could use something like this to reduce copying,
// returning pointers to static data:

  template<T>
  class array_view
  {
  public:
    array_view(std::size_t size, const T * data):
      size(size), data(data)
    const std::size_t size;
    const T * data;
    const T operator[](std::size_t i) const
    { return data[i]; }
  };

  array_view<int> form::original_coefficient_positions() const
  {
    static const int data = { 0, 1, 2, 5 };
    return array_view<int>{4, data};
  }

  array_view<bool> integral::enabled_coefficients() const
  {
    static const bool data = { true, true, true, false, false, true };
    return array_view<bool>{6, data};
  }
"""
