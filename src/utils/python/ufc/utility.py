
def make_unique(li):
    li.sort()
    i = 0
    while i < len(li)-1:
        if li[i+1] == li[i]:
            li.remove(li[i])
        else:
            i += 1
    return li       


def get_format_vars(string):
    vars = []
    pos = string.find("%(")
    while pos >= 0:
        vars.append( string[pos+2:string.find(")s", pos+2)] )
        pos = string.find("%(", pos+1)
    return make_unique(vars)


def indent(string, num_spaces):
    return "\n".join( [" "*num_spaces + li for li in string.split("\n")] )


def generate_code(format_string, format_dict, num_indent_spaces=6):
    """Apply default empty strings for unnecessary variables and indent contents of format_dict."""
    default_dict = {'members':'', 'constructor':'', 'destructor':''}
    default_dict.update(format_dict)
    undefined = []
    vars = get_format_vars(format_string)
    for v in vars:
        default_dict[v] = indent(default_dict[v], num_indent_spaces)
    undefined = make_unique(undefined)
    return format_string % default_dict


