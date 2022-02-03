# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tools for indentation-aware code string stitching.

When formatting an AST into a string, it's better to collect lists of
snippets and then join them than adding the pieces continually, which
gives O(n^2) behaviour w.r.t. AST size n.

"""


class Indented(object):
    """Class to mark a collection of snippets for indentation.

    This way nested indentations can be handled by adding the prefix
    spaces only once to each line instead of splitting and indenting
    substrings repeatedly.

    """

    # Try to keep memory overhead low:
    __slots__ = ("body", )

    def __init__(self, body):
        # Body can be any valid snippet format
        self.body = body


def iter_indented_lines(snippets, level=0):
    """Iterate over indented string lines from a snippets data structure.

    The snippets object can be built recursively using the following types:

    - str: Split and yield as one line at a time indented to the appropriate level.

    - Indented: Yield the lines within this object indented by one level.

    - tuple,list: Yield lines from recursive application of this function to list items.

    """
    tabsize = 2
    indentation = ' ' * (tabsize * level)
    if isinstance(snippets, str):
        for line in snippets.split("\n"):
            yield indentation + line
    elif isinstance(snippets, Indented):
        for line in iter_indented_lines(snippets.body, level + 1):
            yield line
    elif isinstance(snippets, (tuple, list)):
        for part in snippets:
            for line in iter_indented_lines(part, level):
                yield line
    else:
        raise RuntimeError("Unexpected type %s:\n%s" % (type(snippets), str(snippets)))


def format_indented_lines(snippets, level=0):
    """Format recursive sequences of indented lines as one string."""
    return "\n".join(iter_indented_lines(snippets, level))
