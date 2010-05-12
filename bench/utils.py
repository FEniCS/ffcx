__author__ = "Anders Logg"
__date__ = "2010-05-11"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

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

    print ""
    for i in range(m):
        print " " + "-"*row_size
        print "|",
        for j in range(n):
            print table[i][j] + " "*(column_sizes[j] - len(table[i][j])),
            print "|",
        print ""
    print " " + "-"*row_size
    print ""
