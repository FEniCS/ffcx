from ufl.common import product
from six.moves import xrange


class Range(object):

    def __init__(self, begin, end):
        self.begin = begin
        self.end = end


class IndexMapping(object):

    def __init__(self, ranges):
        self.names = sorted(ranges.keys())

        self.ranges = {}
        self.begin = {}
        self.end = {}
        self.size = {}

        for name in self.names:
            r = ranges[name]

            if isinstance(r, int):
                begin, end = 0, r
            elif isinstance(r, str):
                begin, end = "0", r
            elif isinstance(r, Range):
                begin, end = r.begin, r.end
            else:
                begin, end = r

            if begin in (0, "0"):
                size = end
            elif isinstance(end, str) or isinstance(begin, str):
                size = "({0} - {1})".format(end, begin)
            else:
                size = end - begin

            self.ranges[name] = (begin, end)
            self.begin[name] = begin
            self.end[name] = end
            self.size[name] = size

    def __str__(self):
        l = ['{!r}: ({!r}, {!r})'.format(name, self.begin[name], self.end[name]) for name in self.names]
        return '{{ {} }}'.format(', '.join(l))

    def __repr__(self):
        return 'IndexMapping({!s})'.format(self)


def as_tuple(a):
    if isinstance(a, tuple):
        return a
    elif isinstance(a, list):
        return tuple(a)
    else:
        return (a,)


def any_str(*args):
    return any(isinstance(a, str) for a in args)


def dim_mul(i, j):
    if any_str(i, j):
        return '{!r} * {!r}'.format(i, j)
    else:
        return i * j


def mul_dims(dims):
    if not dims:
        return 1
    ints = [d for d in dims if isinstance(d, int)]
    strs = [d for d in dims if isinstance(d, str) and d != "1"]
    if ints and ints != [1]:
        strs = ['{!r}'.format(product(ints))] + strs
    if strs:
        return ' * '.join(strs)
    return '1'


class AxisMapping(object):

    def __init__(self, index_mapping, axes):
        self.index_mapping = index_mapping
        self.axis_index_names = [as_tuple(a) for a in axes]

        self.num_axes = len(self.axis_index_names)
        self.num_dims = [len(a) for a in self.axis_index_names]

        # The size of each axis equals the product of the associated index dimensions
        self.axis_size = [mul_dims([self.index_mapping.size[name] for name in names])
                          for names in self.axis_index_names]

        # The stride of each axis equals the product of the following axis sizes
        self.axis_stride = [mul_dims(self.axis_size[i + 1:]) for i in range(self.num_axes)]

        # For each axis, the internal strides of each dimension within that axis
        self.dim_stride = [tuple(mul_dims([self.index_mapping.size[name] for name in self.axis_index_names[i][j + 1:]])
                                 for j in range(self.num_dims[i]))
                           for i in range(self.num_axes)]

    def format_decl(self):
        return ''.join('[{0}]'.format(size) for size in self.axis_size)

    def format_access(self, **kwargs):
        expressions = self.format_index_expressions(**kwargs)
        return ''.join('[{0}]'.format(expr) for expr in expressions)

    def format_index_expressions(self, **kwargs):
        expressions = []
        for i in range(self.num_axes):
            terms = []
            for j in range(self.num_dims[i]):
                stride = self.dim_stride[i][j]
                index = kwargs.get(self.axis_index_names[i][j], self.axis_index_names[i][j])
                term = mul_dims([stride, index])
                terms.append(term)
            expr = ' + '.join(terms)
            expressions.append(expr)
        return expressions

# TODO: Change this into unit tests


def _test():
    ranges = {"i": 4, "j": (10, 12), "n": "N", "m": (1, "M")}
    im = IndexMapping(ranges)
    print(str(eval(repr(im))))

    axes = ["i", "j", "n", "m"]
    am = AxisMapping(im, axes)
    print()
    print(am.axis_index_names)
    print(am.axis_size)
    print(am.axis_stride)
    print(am.dim_stride)
    print(am.format_decl())
    print(am.format_access(i='ii', j='jj', n='nn', m='mm'))

    axes = [("i", "j", "n", "m")]
    am = AxisMapping(im, axes)
    print()
    print(am.axis_index_names)
    print(am.axis_size)
    print(am.axis_stride)
    print(am.dim_stride)
    print(am.format_decl())
    print(am.format_access(i='ii', j='jj', n='nn', m='mm'))

#_test()
