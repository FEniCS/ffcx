#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

class test_codegen_testing_framework(CodegenTestCase):
    '''The part of this docstring above the capital header: line
    will not be included in the generated test code. The part
    below the capital header: line will be inserted at the top
    of the test code.

    Docstrings for each test function will be split into three parts.
    The initial part becomes the docstring for the C++ test code.
    The part after the PRE: line is inserted before the code under test.
    The part after the POST: line is inserted after the code under test.

    HEADER:
    /**
    Unit tests of CodegenTestCase framework.
    */
    #include <iostream>
    using std::cout;
    using std::endl;
    '''

    def test_empty_docstring(self):
        """
        """
        code = '    ASSERT_TRUE(1);'
        # Enable this line to see screen output from the generated test:
        #code += '\n    cout << "Hello Test World!" << endl;'
        self.emit_test(code)

    def test_nopre_nopost_docstring(self):
        """Documentation.
        """
        code = "    ASSERT_TRUE(1);"
        self.emit_test(code)

    def test_nopre_docstring(self):
        """Documentation.

        POST:
        state = true;  // post code 1/2
        ASSERT_TRUE(state); // post code 2/2
        """
        code = "    bool state = false;"
        self.emit_test(code)

    def test_nopost_docstring(self):
        """Documentation.

        PRE:
        bool state = false;  // pre code 1/2
        state = true;        // pre code 2/2
        """
        code = "    ASSERT_TRUE(state);"
        self.emit_test(code)

    def test_full_docstring(self):
        """Documentation 1/2.
        Documentation 2/2.

        PRE:
        bool state = false;  // pre code 1/2
        ; // pre code 2/2

        POST:
        ASSERT_TRUE(state); // post code 1/2
        ; // post code 2/2
        """
        code = "    state = true;"
        self.emit_test(code)

    def test_compact_full_docstring(self):
        """PRE:
        POST:
        """
        code = "    ASSERT_TRUE(1);"
        self.emit_test(code)

    def test_compact_full_docstring2(self):
        """Doc.
        PRE:
        bool state = false; // pre code
        POST:
        ASSERT_TRUE(state); // post code 1/1
        """
        code = "    state = true;"
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
