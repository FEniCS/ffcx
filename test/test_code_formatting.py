#!/usr/bin/env python

"""
"""

# These are thin wrappers on top of unittest.TestCase and unittest.main
from ufltestcase import UflTestCase, main

class CodeFormattingTestCase(UflTestCase):

    def test_c_code_formatting(self):
        from uflacs.c_format_test import test_code_formatting
        test_code_formatting()

    def test_format_code(self):
        from uflacs.format_code import test_format_code
        test_format_code()


if __name__ == "__main__":
    main()

