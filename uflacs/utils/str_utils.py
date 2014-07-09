"""
String manipulation and formatting utilities.
"""

from six.moves import iteritems

def format_list(obj):
    assert isinstance(obj, (list, tuple))
    return '\n'.join('  %s' % str(v) for v in obj)

def format_dict(obj):
    assert isinstance(obj, dict)
    fmt = lambda k, v: '  %s: %s' % (k, v)
    it = sorted(iteritems(obj), key=lambda x: x[0])
    return '\n'.join(['{', '\n'.join(fmt(k, v) for (k,v) in it), '}'])

def format_obj(obj):
    if isinstance(obj, (list, tuple)):
        return format_list(obj)
    if isinstance(obj, dict):
        return format_dict(obj)
    return str(obj)

def format_members(obj, prefix=''):
    fmt = lambda k, v: '%s.%s = \n%s' % (prefix, k, format_obj(v))
    it = sorted(iteritems(vars(obj)), key=lambda x: x[0])
    return "\n\n".join(fmt(k, v) for (k,v) in it)


# TODO: These are basically the same as the above:

def format_sequence(sequence):
    return '\n'.join("{0}".format(v) for v in sequence)

def format_enumerated_sequence(sequence):
    return '\n'.join("{0}: {1}".format(i, v) for i,v in enumerate(sequence))

def format_mapping(mapping):
    return '\n'.join("{0}: {1}".format(k, v) for k,v in mapping.items())
