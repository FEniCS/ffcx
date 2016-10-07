#!/usr/bin/env py.test
# -*- coding: utf-8 -*-

def test_example_showing_how_to_test_generated_code_with_gtest(gtest):
    "This is an example test explaining the py.test/gtest integration framework."

    # First we 'generate' some code
    cut = "x[0] = 0.0; x[1] = 2.0;"

    # We can make regular py.test checks as we go
    assert "x[0] = 0.0;" in cut

    # Then we wrap the generated code in C++ setup and gtest checks
    gtestcode = """
    double x[2] = {{ 1.0, 2.0 }};
    {cut}
    ASSERT_TRUE(x[0] == 0.0);
    """.format(cut=cut)

    # Finally we emit the C++ test, which will be run in a later pass
    gtest.add(gtestcode)
