# -*- coding: utf-8 -*-

from ffc.uflacs.language.format_value import format_float
from ffc.uflacs.language.format_lines import iter_indented_lines, Indented, format_indented_lines


def test_format_float():
    # Ints handled
    assert format_float(0, 3) == "0.0"
    assert format_float(1, 3) == "1.0"
    assert format_float(12, 15) == "12.0"

    # Zeros simple
    assert format_float(0.0, 0) == "0.0"
    assert format_float(0.0, 3) == "0.0"
    assert format_float(0.0, 15) == "0.0"

    # Ones simple
    assert format_float(1.0, 0) == "1.0"
    assert format_float(1.0, 3) == "1.0"
    assert format_float(1.0, 15) == "1.0"

    # Small ints simple
    assert format_float(12., 15) == "12.0"
    assert format_float(12., 15) == "12.0"

    # Precision truncates
    assert format_float(1.2, 3) == "1.2"
    assert format_float(1.23, 3) == "1.23"
    assert format_float(1.2345, 3) == "1.23"
    assert format_float(1.0 + 1e-5, 7) == "1.00001"
    assert format_float(1.0 + 1e-5, 6) == "1.00001"
    assert format_float(1.0 + 1e-5, 5) == "1.0"

    # Cleanly formatted exponential numbers
    # without superfluous +s and 0s
    assert format_float(1234567.0, 3) == "1.23e6"
    assert format_float(1.23e6, 3) == "1.23e6"
    assert format_float(1.23e-6, 3) == "1.23e-6"
    assert format_float(-1.23e-6, 3) == "-1.23e-6"
    assert format_float(-1.23e6, 3) == "-1.23e6"



def test_iter_indented_lines():
    assert list(iter_indented_lines("single line")) == ["single line"]
    assert list(iter_indented_lines(["word"])) == ["word"]
    assert list(iter_indented_lines(["line one", "line two"])) == ["line one", "line two"]
    assert list(iter_indented_lines("line one\nline two")) == ["line one", "line two"]
    assert list(iter_indented_lines(["line one", Indented("line two")])) == ["line one", "    line two"]
    assert list(iter_indented_lines([Indented("line one"), "line two"])) == ["    line one", "line two"]
    assert list(iter_indented_lines([[Indented("line one")], "line two"])) == ["    line one", "line two"]
    assert list(iter_indented_lines([Indented("line one"), ["line two"]])) == ["    line one", "line two"]
    assert list(iter_indented_lines([[Indented("line one"), "line two"]])) == ["    line one", "line two"]


def test_format_indented_lines_example():
    forheader = "for (int i=begin; i!=end; ++i)"
    code = ["{",
            Indented([forheader,
                      "{",
                      Indented("double x = y[i];"),
                      "}"]),
            "}"]
    reference_code = """{
    for (int i=begin; i!=end; ++i)
    {
        double x = y[i];
    }
}"""
    assert format_indented_lines(code) == reference_code
