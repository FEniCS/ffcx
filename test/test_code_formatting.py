#!/usr/bin/env python

"""
"""

# These are thin wrappers on top of unittest.TestCase and unittest.main
from ufltestcase import UflTestCase, main

class CodeFormattingTestCase(UflTestCase):

    def test_foo(self):
        from uflacs.c_format_test import test_code_formatting
        test_code_formatting()

if __name__ == "__main__":
    main()

