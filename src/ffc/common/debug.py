"""A common interface for diagnostic messages. Only messages with
debug level lower than or equal to the current debug level will be
printed. To see more messages, raise the debug level."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-02-04 -- 2007-02-06"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

import sys

__level = -1
__indent = 0
__no_indent = False

def debug(string, debuglevel = 0):
    "Print given string at given debug level"

    global __level, __no_indent

    # Print message (not fancy handling of strings containing ...)
    if debuglevel <= __level:
        indentation = "".join(["    " for i in range(__indent)])
        if "..." in string:
            print indentation + string,
            __no_indent = True
        elif __no_indent:
            print string
            __no_indent = False
        else:
            print indentation + string

    # Flush buffer so messages are printed *now*
    sys.stdout.flush()

def debug_indent(increment = 1):
    "Set indentation of debug messages"
    global __indent
    __indent += increment

def setlevel(newlevel):
    "Set debug level"
    global __level
    __level = newlevel

def getlevel():
    "Get debug level"
    return __level
