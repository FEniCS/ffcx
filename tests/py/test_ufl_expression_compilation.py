#!/usr/bin/env python
from codegentestcase import CodegenTestCase, unittest

import ufl
from ufl.common import product

from uflacs.backends.cpp2.compiler import compile_expression

class test_ufl_expression_compilation(CodegenTestCase):
    '''TODO: Fill in the blanks!

    HEADER:
    /**
    Unit tests of generated geometry snippet code.
    */
    #include <iostream>
    using std::cout;
    using std::endl;
    '''

    def compile_expression0(self, expr):
        code = ""
        return code

    def compile_expression1(self, expr):
        code = "double values[%d];" % (product(expr.shape()),)
        return code

    def compile_expression(self, expr):
        code = compile_expression(expr, "")
        return code

    def test_compilation_of_x(self):
        """Test that compilation of x results in the value of x placed in values.

        PRE:
        double x[3] = { 0.0, 0.1, 0.2 };
        double A[3];

        POST:
        ASSERT_EQ(A[0], 0.0);
        ASSERT_EQ(A[1], 0.1);
        ASSERT_EQ(A[2], 0.2);
        """
        expr = ufl.tetrahedron.x
        code = self.compile_expression(expr)
        self.emit_test(code)

if __name__ == "__main__":
    unittest.main()
