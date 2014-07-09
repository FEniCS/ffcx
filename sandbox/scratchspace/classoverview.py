
from six import iteritems
from six.moves import map

def subclass_tree(base, all_classes):
    # This is O(n**2), but I don't care
    subtrees = {}
    for c in all_classes:
        if c.__base__ is base:
            subtrees[c] = subclass_tree(c, all_classes)
    return subtrees

def format_tree_indented(tree, indent=0, formatter=str):
    code = []
    ind = '\t'*indent
    def cmpformatted(x, y):
        return cmp(formatter(x), formatter(y))
    for k in sorted(tree, cmp=cmpformatted):
        v = tree[k]
        lineend = ':' if v else ''
        code.append('%s%s%s' % (ind, formatter(k), lineend))
        if v:
            code.append(format_tree_indented(v, indent=indent+1, formatter=formatter))
    return '\n'.join(c for c in code)

def flatten_tree(tree):
    flat = []
    for node in sorted(tree):
        flat.append(node)
        flat.extend(flatten_tree(tree[node]))
    return flat

def get_callables(module):
    name = module.__name__
    def cond(v):
        return callable(v) and v.__module__.startswith(name)
    return [n for n, v in iteritems(vars(module)) if cond(v)]

def curry(f, g):
    def h(*args, **kwargs):
        return f(g(*args, **kwargs))
    return h

def getname(x):
    return x.__name__

def indent(x):
    return '\t' + x

def main():
    from ufl.classes import Expr, Terminal, Operator, all_ufl_classes, abstract_classes, terminal_classes
    tree = subclass_tree(Expr, all_ufl_classes)
    print(format_tree_indented(tree, formatter=getname))

    # Find classes most interesting to test:
    flat = flatten_tree(tree)
    nonabstract = [c for c in flat if c not in abstract_classes]
    terminal    = [c for c in nonabstract if issubclass(c, Terminal)]
    nonterminal = [c for c in nonabstract if not isinstance(c, Terminal)]

    # Print classes
    fmt = curry(indent, getname)
    fmt = getname
    print("\nNonabstract UFL terminal types:")
    print('\n'.join(map(fmt, terminal)))
    print("\nNonabstract UFL operator types:")
    print('\n'.join(map(fmt, nonterminal)))
    print('\n%d terminal types, %d operator types' % (len(terminal), len(nonterminal)))

    # Print algorithms
    import ufl.algorithms
    callables = get_callables(ufl.algorithms)
    print('\n'.join(sorted(callables)))
    print('\n%d functions from algorithms' % len(callables))

if __name__ == '__main__':
    main()

# Reminder to self:
# TODO: min/max in ufl?
# TODO: atan2 in ufl?
