
# TODO: Move to uflacs.language utils?
def generate_return_new_switch(L, i, classnames):
    if classnames:
        cases = []
        for j, classname in enumerate(classnames):
            if classname:
                cases.append((j, L.Return(L.New(classname))))
        code = [L.Switch(i, cases, autobreak=False, autoscope=False)]
    else:
        code = []
    code.append(L.Return(L.Null()))
    return L.StatementList(code)
