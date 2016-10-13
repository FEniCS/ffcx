# -*- coding: utf-8 -*-

from ffc.uflacs.language.format_value import format_float, set_float_precision, reset_float_precision
from ffc.uflacs.language.format_lines import iter_indented_lines, Indented, format_indented_lines


def test_format_float():
    reset_float_precision()
    assert format_float(0.0) == "0.0"
    assert format_float(1.0) == "1.0"
    assert format_float(12.) == "12.0"

    set_float_precision(3)
    assert format_float(0.0) == "0.0"
    assert format_float(1.0) == "1.0"
    assert format_float(1.2) == "1.2"
    assert format_float(1.23) == "1.23"

    set_float_precision(15, 1e-15)
    assert format_float(0.0) == "0.0"
    assert format_float(1.0) == "1.0"
    assert format_float(12.) == "12.0"  # 1.2e+01

    reset_float_precision()


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
