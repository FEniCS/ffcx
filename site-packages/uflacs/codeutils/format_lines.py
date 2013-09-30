
def format_lines(fmt, lineargs):
    for args in lineargs:
        yield fmt % args

def format_assignments(name_value_pairs):
    fmt = "%s = %s;"
    return format_lines(fmt, name_value_pairs)

def format_additions(name_value_pairs):
    fmt = "%s += %s;"
    return format_lines(fmt, name_value_pairs)

def format_scaled_assignments(name_value_pairs, factor):
    fmt = "%%s = %s * %%s;" % factor
    return format_lines(fmt, name_value_pairs)

def format_scaled_additions(name_value_pairs, factor):
    fmt = "%%s += %s * %%s;" % factor
    return format_lines(fmt, name_value_pairs)
