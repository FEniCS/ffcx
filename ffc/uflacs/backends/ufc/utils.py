# -*- coding: utf-8 -*-


# TODO: Move these to uflacs.language utils?


def generate_return_new(L, classname, factory):
    if factory:
        return L.Return(L.Call("create_" + classname))
    else:
        return L.Return(L.New(classname))


def generate_return_new_switch(L, i, classnames, args=None, factory=False):
    if factory:
        def create(classname):
            return L.New(classname)
    else:
        def create(classname):
            return L.Call("create_" + classname)
    if classnames:
        cases = []
        if args is None:
            args = list(range(len(classnames)))
        for j, classname in zip(args, classnames):
            if classname:
                cases.append((j, L.Return(create(classname))))
        code = [L.Switch(i, cases, autobreak=False, autoscope=False)]
    else:
        code = []
    code.append(L.Return(L.Null()))
    return L.StatementList(code)
