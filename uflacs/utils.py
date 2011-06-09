
def format_list(obj):
    assert isinstance(obj, (list, tuple))
    return '\n'.join('  %s' % str(v) for v in obj)

def format_dict(obj):
    assert isinstance(obj, dict)
    format = lambda k, v: '  %s = %s' % (k, v)
    it = obj.iteritems()
    return '\n'.join(format(k, v) for (k,v) in it)

def format_obj(obj):
    if isinstance(obj, (list, tuple)):
        return format_list(obj)
    if isinstance(obj, dict):
        return format_dict(obj)
    return str(obj)

def format_members(obj):
    format = lambda k, v: '%s:\n%s' % (k, format_obj(v))
    it = vars(obj).iteritems()
    return "\n\n".join(format(k,v) for (k,v) in it)
