# -*- coding: utf-8 -*-

from six import string_types

# TODO: Move these to uflacs.language utils?


def generate_return_new(L, classname, factory):
    if factory:
        return L.Return(L.Call("create_" + classname))
    else:
        return L.Return(L.New(classname))


def generate_return_new_switch(L, i, classnames, args=None, factory=False):
    # TODO: UFC functions of this type could be replaced with return vector<shared_ptr<T>>{objects}.

    if isinstance(i, string_types):
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


def generate_return_literal_switch(L, i, values, default, literal_type):
    # TODO: UFC functions of this type could be replaced with return vector<T>{values}.

    if isinstance(i, string_types):
        i = L.Symbol(i)

    default = L.Return(literal_type(default))
    if values:
        cases = [(j, L.Return(literal_type(k)))
                 for j, k in enumerate(values)]
        return L.Switch(i, cases, default=default)
    else:
        return default


def generate_return_int_switch(L, i, values, default):
    return generate_return_literal_switch(L, i, values, default, int)


def generate_return_bool_switch(L, i, values, default):
    return generate_return_literal_switch(L, i, values, default, bool)


# TODO: Better error handling
def generate_error(L, msg, emit_warning):
    if emit_warning:
        return L.VerbatimStatement('std::cerr << "*** FFC warning: " << "%s" << std::endl;' % (msg,))
    else:
        return L.Raise('std::runtime_error("%s");' % (msg,))
