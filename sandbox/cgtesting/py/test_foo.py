#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

class test_foo(CodegenTestCase):
    '''This is just a dummy test.

    HEADER:
    /**
    Unit tests of CodegenTestCase framework.
    */
    #include <iostream>
    using std::cout;
    using std::endl;
    '''

    def test_some(self):
        code = 'double some() { cout << "some" << endl; return 0.1; }'
        self.emit_code(code)

    def test_some_more(self):
        code = 'double some_more() { cout << "some more" << endl; return 0.2; }'
        self.emit_code(code)

if __name__ == "__main__":
    unittest.main()
