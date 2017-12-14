# -*- coding: utf-8 -*-
# Copyright (C) 2010 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2010-05-11
# Last changed: 2010-05-11

def print_table(values, title):
    "Print nicely formatted table."

    m = max([key[0] for key in values]) + 2
    n = max([key[1] for key in values]) + 2

    table = []
    for i in range(m):
        table.append(["" for j in range(n)])

    for i in range(m - 1):
        table[i + 1][0] = str(values[(i, 0)][0])

    for j in range(n - 1):
        table[0][j + 1] = str(values[(0, j)][1])

    for i in range(m - 1):
        for j in range(n - 1):
            value = values[(i, j)][2]
            if isinstance(value, float):
                value = "%.5g" % value
            table[i + 1][j + 1] = value

    table[0][0] = title

    column_sizes = [max([len(table[i][j]) for i in range(m)]) for j in range(n)]
    row_size = sum(column_sizes) + 3*(len(column_sizes) - 1) + 2

    print("")
    for i in range(m):
        print(" " + "-"*row_size)
        print("|", end="")
        for j in range(n):
            print(table[i][j] + " "*(column_sizes[j] - len(table[i][j])), end="")
            print("|", end="")
        print("")
    print(" " + "-"*row_size)
    print("")
