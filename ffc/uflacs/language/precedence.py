# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

class PRECEDENCE:
    "An enum-like class for C operator precedence levels."
    HIGHEST = 0
    LITERAL = 0
    SYMBOL = 0

    #SCOPE = 1

    POST_INC = 2
    POST_DEC = 2
    CALL = 2
    SUBSCRIPT = 2
    #MEMBER = 2
    #PTR_MEMBER = 2

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

    BIT_AND = 8
    BIT_XOR = 9
    BIT_OR = 10
    AND = 11
    OR = 12

    CONDITIONAL = 13
    ASSIGN = 13

    #COMMA = 14

    LOWEST = 15
