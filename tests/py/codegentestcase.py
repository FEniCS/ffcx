import unittest
import os
from glob import glob
import inspect
from itertools import takewhile

main_template = """
/*******************************
 *
 *  Generated Google test framework main file.
 *
 ********************************/

#include <gtest/gtest.h>

%s

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
"""

gc_header = """
/*******************************
 *
 *  Generated code ready for testing.
 *  NB! This file is regenerated each
 *  run, so manual edits will be lost.
 *
 ********************************/
"""

test_skeleton_template = """\
#ifndef %(upper)s_H_INCLUDED
#define %(upper)s_H_INCLUDED

#include <gtest/gtest.h>

#include "%(gcname)s"

TEST (%(basename)s, test_dummy)
{
    ASSERT_EQ(1, 1);
}

#endif
"""

test_case_template = """/**
%s
*/
TEST (%s, %s)
{
    // Precondition code:
%s

    // Generated code:
%s

    // Postcondition code:
%s
}"""

def indented(lines):
    indent = '    '
    return '\n'.join(indent+l for l in lines)

def format_test_case(doc, suite, case, pre, code, post):
    code = code.rstrip().lstrip('\n')
    return test_case_template % (doc, suite, case, pre, code, post)

def format_includes(includes, system_includes):
    return '\n'.join(['\n'.join('#include "%s"' % inc for inc in includes),
                      '\n'.join('#include <%s>' % inc for inc in system_includes)])

def format_test_skeleton(basename, gcname):
    upper = basename.upper()
    return test_skeleton_template % locals()

def format_main(includes):
    return main_template % format_includes(includes, [])

def split_by_line_markers(text, markers):
    lines = [l.strip('\n') for l in text.split('\n')]

    # Run through all lines and keep the current
    # paragraph name in active_marker, adding
    # named paragraphs at each change of active_marker
    active_marker = ""
    paragraphs = {}
    paragraph = []
    for line in lines:
        m = line.strip()
        if m in markers:
            assert active_marker not in paragraphs
            paragraphs[active_marker] = indented(paragraph).rstrip()

            paragraph = []
            active_marker = m
        else:
            paragraph.append(line)

    assert active_marker not in paragraphs
    paragraphs[active_marker] = indented(paragraph).rstrip()

    return [paragraphs.get(m,"") for m in [""]+list(markers)]

def find_parent_test_function():
    frame = inspect.currentframe()
    function = ""
    while not function.startswith("test_"):
        frame = frame.f_back
        function = inspect.getframeinfo(frame)[2]
    return function

def basename2gcname(base_name):
    return "gc_" + base_name + ".h"

def create_initial_test_header(tfn, bn=None, ghn=None):
    # Extract names from test header name, because this
    # can be run from commandline as well as from test class
    base_name = os.path.basename(tfn).replace(".h", "")
    gc_header_name = basename2gcname(base_name)
    assert bn is None or bn == base_name
    assert ghn is None or ghn == gc_header_name

    f = open(tfn, 'w')
    f.write(format_test_skeleton(base_name, gc_header_name))
    f.close()

class CodegenTestCase(unittest.TestCase):
    visited = set()

    # NB! These paths should match the makefile in the tests directory
    gendir = "generated"
    cppdir = "cpp"

    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)

        # Set to nonzero to show info about created files and emitted code
        self.verbose = 1

        # Names specific to each subclass
        self.base_name = type(self).__name__
        self.test_header_name = self.base_name + ".h"
        self.gc_header_name = basename2gcname(self.base_name)
        self.main_src_name = "main_" + self.base_name + ".cpp"

        self.visit()

    def visit(self):
        if self.base_name not in CodegenTestCase.visited:
            # Keep a record for this python session, remembering
            # which files have been initialized. Since each emit will
            # be run in a separate instance of a test suite class,
            # the visited set must be stored as a "singleton".
            CodegenTestCase.visited.add(self.base_name)

            # Generate top part of header file for emitted code,
            # this will be regenerated each run
            gfn = os.path.join(self.gendir, self.gc_header_name)
            f = open(gfn, 'w')
            f.write(gc_header)
            doc, header = split_by_line_markers(self.__doc__, ("HEADER:",))
            f.write(header)
            f.close()
            if self.verbose:
                print "Initialized", gfn

            # Generate initial test header, later this can
            # be manually modified so don't overwrite it
            tfn = os.path.join(self.cppdir, self.test_header_name)
            if not glob(tfn):
                create_initial_test_header(tfn,
                                           self.base_name, self.gc_header_name)
                if self.verbose:
                    print "Initialized", tfn

            # Generate initial main src file for this suite, later
            # this can be manually modified so don't overwrite it
            mfn = os.path.join(self.gendir, self.main_src_name)
            if not glob(mfn):
                f = open(mfn, 'w')
                f.write(format_main([self.test_header_name]))
                f.close()
                if self.verbose:
                    print "Initialized", mfn

    def emit_test(self, code):
        # Look through stack to fine frame with test_... function name
        function = find_parent_test_function()

        # Using testcase class name as suite and function name as case
        suite = self.__class__.__name__
        case = function

        # Extracting C++ testcase docstring, precondition code, and
        # postcondition code from docstring of python test function
        fulldoc = inspect.getdoc(getattr(self, function)) or "Generated test with docstring."
        doc, pre, post = split_by_line_markers(fulldoc, ("PRE:", "POST:"))

        # Combine snippets into a full test case
        testcode = format_test_case(doc, suite, case, pre, code, post)

        if self.verbose:
            print "EMITTING TEST ", suite, case

        # Output test case code to file included in C++ tests
        self.emit_code(testcode)

    def emit_code(self, code):
        # Output code directly to file included in C++ tests
        hfn = os.path.join(self.gendir, self.gc_header_name)
        f = open(hfn, 'a')
        f.write('\n')
        f.write(code)
        f.write('\n')
        f.close()
        if self.verbose:
            print "Wrote code to file", hfn
