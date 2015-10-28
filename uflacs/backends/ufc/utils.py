
# TODO: Move to uflacs.language utils?
def generate_return_new_switch(L, i, classnames, args=None):
    if classnames:
        cases = []
        if args is None:
            args = list(range(len(classnames)))
        for j, classname in zip(args, classnames):
            if classname:
                cases.append((j, L.Return(L.New(classname))))
        code = [L.Switch(i, cases, autobreak=False, autoscope=False)]
    else:
        code = []
    code.append(L.Return(L.Null()))
    return L.StatementList(code)
