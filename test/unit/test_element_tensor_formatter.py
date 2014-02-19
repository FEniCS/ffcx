#!/usr/bin/env python
"""
Tests of generic C++ compilation code.
"""

import uflacs
from uflacs.codeutils.format_lines import *
from uflacs.codeutils.element_tensor_formatter import build_loops

import ufl
from ufl import *
from uflacs.codeutils.format_code_structure import format_code_structure

def test_format_lines():
    def lineargs1():
        return ()
    fmt = "dummy"
    r = list(format_lines(fmt, lineargs1()))
    assert r == []

    def lineargs2():
        for i in range(2):
            yield "hello", i*2, i*3
    fmt = "%s, %d, %d"
    r = list(format_lines(fmt, lineargs2()))
    assert r == ["hello, 0, 0", "hello, 2, 3"]

def test_format_assignments():
    def name_value_pairs():
        for i in range(2):
            yield "s[%d]" % i, str((i+1)**2)
    r = list(format_assignments(name_value_pairs()))
    assert r == ["s[0] = 1;", "s[1] = 4;"]

def test_format_additions():
    def name_value_pairs():
        for i in range(2):
            yield "A[%d]" % i, str((i+1)**2)
    r = list(format_additions(name_value_pairs()))
    assert r == ["A[0] += 1;", "A[1] += 4;"]

def test_build_loops_empty():
    loops = []
    definitions = []
    partitions = []
    code = build_loops(loops, definitions, partitions)
    expected = ""
    assert format_code_structure(code) == expected

def test_build_loops_one_level():
    """Build a single loop."""
    loops = ["for (i0;;)"]
    definitions = ["defs0;"]
    partitions = ["part0;"]
    code = build_loops(loops, definitions, partitions)
    #print format_code_structure(code)
    expected = "for (i0;;)\n{\n    defs0;\n    part0;\n}"
    assert format_code_structure(code) == expected

def test_build_loops_two_levels_two_loops():
    """Build nested code with inner and outer loop."""
    loops = ["for (i0;;)", "for (i1;;)"]
    definitions = ["defs0;", "defs1;"]
    partitions = ["part0;", "part1;"]
    code = build_loops(loops, definitions, partitions)
    #print format_code_structure(code)
    expected = "for (i0;;)\n{\n    defs0;\n    part0;\n    "+\
        "for (i1;;)\n    {\n        defs1;\n        part1;\n    }\n}"
    assert format_code_structure(code) == expected

def test_build_loops_two_levels_only_inner_loop():
    """Build nested code with no outer loop."""
    loops = ["", "for (i1;;)"]
    definitions = ["defs0;", "defs1;"]
    partitions = ["part0;", "part1;"]
    code = build_loops(loops, definitions, partitions)
    #print format_code_structure(code)
    expected = "defs0;\npart0;\nfor (i1;;)\n{\n    defs1;\n    part1;\n}"
    assert format_code_structure(code) == expected

def test_build_loops_two_equivalent_levels():
    """If you want to make two loops after each other,
    just use build_loops twice."""
    code = []

    loops = ["for (i0;;)"]
    definitions = ["defs0;"]
    partitions = ["part0;"]
    code.append(build_loops(loops, definitions, partitions))

    loops = ["for (i1;;)"]
    definitions = ["defs1;"]
    partitions = ["part1;"]
    code.append(build_loops(loops, definitions, partitions))

    #print format_code_structure(code)
    expected = "for (i0;;)\n{\n    defs0;\n    part0;\n}\nfor (i1;;)\n{\n    defs1;\n    part1;\n}"
    assert format_code_structure(code) == expected

def test_toy_compiler():
    # This code could become the basis of a backend-generic form compiler

    # Assignments simulating partitioned computational graph
    p0_assignments = [("p0[0]", "3.14")]
    p1_assignments = [("p1[0]", "p0[0]*x[0]")]
    p2_assignments = [("p2[0]", "p1[0]*v0")]
    A_accumulations = [("A[0]", "p2[0]*v1")]

    defs_cell = ["int nq = 1;",
                 "int n0 = 1;",
                 "int n1 = 1;"]
    part_cell = ["double p0[1];",
                 list(format_assignments(p0_assignments))]

    loop_quad = "for (int iq=0; iq<nq; ++iq)"
    defs_quad = ["// Define geometry in x[iq]",
                 "double x = xq[iq];"]
    part_quad = ["double p1[1];",
                 list(format_assignments(p1_assignments))]

    loop_test = "for (int i0=0; i0<n0; ++i0)"
    defs_test = ["// Evaluate test function i0 at x[iq]",
                 "double v0 = v0q[i0][i0];"]
    part_test = ["double p2[1];",
                 list(format_assignments(p2_assignments))]

    loop_trial = "for (int i1=0; i1<n1; ++i1)"
    defs_trial = ["// Evaluate trial function i1 at x[iq]",
                  "double v1 = v1q[i1][i0];"]
    part_trial = ["// Accumulate into element tensor",
                  list(format_additions(A_accumulations))]

    loops = ["", loop_quad, loop_test, loop_trial]
    definitions = [defs_cell, defs_quad, defs_test, defs_trial]
    partitions = [part_cell, part_quad, part_test, part_trial]
    code = build_loops(loops, definitions, partitions)

    body = format_code_structure(code)
    #print body
    # Should maybe assert something here, only checked by manual inspection
