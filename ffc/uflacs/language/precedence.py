# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later


class PRECEDENCE:
    """An enum-like class for C operator precedence levels."""

    HIGHEST = 0
    LITERAL = 0
    SYMBOL = 0

    # SCOPE = 1

    POST_INC = 2
    POST_DEC = 2
    CALL = 2
    SUBSCRIPT = 2
    # MEMBER = 2
    # PTR_MEMBER = 2

    PRE_INC = 3
    PRE_DEC = 3
    NOT = 3
    BIT_NOT = 3
    POS = 3
    NEG = 3
    DEREFERENCE = 3
    ADDRESSOF = 3
    SIZEOF = 3

    MUL = 4
    DIV = 4
    MOD = 4

    ADD = 5
    SUB = 5

    LT = 6
    LE = 6
    GT = 6
    GE = 6

    EQ = 7
    NE = 7

    AND = 11
    OR = 12

    CONDITIONAL = 13
    ASSIGN = 13

    # COMMA = 14

    LOWEST = 15
