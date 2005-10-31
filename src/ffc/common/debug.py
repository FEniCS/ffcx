__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-02-04 -- 2005-10-30"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

__level = 0
__indent = 0

"""Diagnostic messages are passed through a common interface to make
it possible to turn debugging on or off. Only messages with debug
level lower than or equal to the current debug level will be
printed. To see more messages, raise the debug level."""

def debug(string, debuglevel = 0):
    global __level
    if debuglevel <= __level:
        indentation = "".join(["    " for i in range(__indent)])
        print indentation + string

def level(newlevel):
    global __level
    __level = newlevel

#def indent(increment = 1):
#    global __indent
#    __indent += increment
