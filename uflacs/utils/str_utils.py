"""
String manipulation and formatting utilities.
"""

def format_list(obj):
    assert isinstance(obj, (list, tuple))
    return '\n'.join('  %s' % str(v) for v in obj)

def format_dict(obj):
    assert isinstance(obj, dict)
    fmt = lambda k, v: '  %s: %s' % (k, v)
    it = sorted(obj.iteritems(), key=lambda x: x[0])
    return '\n'.join(['{', '\n'.join(fmt(k, v) for (k,v) in it), '}'])

def format_obj(obj):
    if isinstance(obj, (list, tuple)):
        return format_list(obj)
    if isinstance(obj, dict):
        return format_dict(obj)
    return str(obj)

def format_members(obj, prefix=''):
    fmt = lambda k, v: '%s.%s = \n%s' % (prefix, k, format_obj(v))
    it = sorted(vars(obj).iteritems(), key=lambda x: x[0])
    return "\n\n".join(fmt(k, v) for (k,v) in it)

