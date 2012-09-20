#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

class MockBasicElementCodeGenerator:
    pass

class test_element_combinations(CodegenTestCase):
    '''TODO: Fill in the blanks!

    HEADER:
    /**
    Unit tests of generated element hierarchy combination code.
    */
    #include <iostream>
    using std::cout;
    using std::endl;
    '''

    def test_mock_element_generator(self):
        """Test that the basic element mock class works as assumed.

        PRE:
        double x[3] = { 0.0, 0.1, 0.2 };

        POST:
        ASSERT_EQ(xi[0], 0.0);
        ASSERT_EQ(xi[1], 0.1);
        ASSERT_EQ(xi[2], 0.1);
        """

        # FIXME: Define interface for mock element code generator using TDD

        code = """
        double xi[3];
        xi[0] = x[0];
        xi[1] = x[1];
        xi[2] = x[2];
        """

        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
