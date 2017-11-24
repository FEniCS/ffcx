# -*- coding: utf-8 -*-
"This script plots the results found in bench.log."

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
# First added:  2010-05-13
# Last changed: 2010-05-13

from pylab import *

# Read logfile
results = {}
try:
    output = open("bench.log").read()
except Exception:
    output = open("results/bench.log").read()
for line in output.split("\n"):
    if "," not in line: continue
    test_case, test_option, timing = [w.strip() for w in line.split(",")]
    try:
        form, degree = test_case.split("_")
    except Exception:
        form, dim, degree = test_case.split("_")
        form = form + "_" + dim
    if form not in results:
        results[form] = {}
    if test_option not in results[form]:
        results[form][test_option] = ([], [])
    results[form][test_option][0].append(int(degree))
    results[form][test_option][1].append(float(timing))

# Plot results
forms = sorted([form for form in results])
test_options = ["-r quadrature", "-r quadrature -O", "-r tensor", "-r tensor -O"]
bullets = ["x-", "o-", "*-", "s-"]
for (i, form) in enumerate(forms):
    figure(i)

    # Plot timings
    subplot(121)
    for (j, test_option) in enumerate(test_options):
        q, t = results[form][test_option]
        semilogy(q, t, bullets[j])
        hold(True)
    a = list(axis()); a[-1] = 50.0*a[-1]; axis(a);
    legend(test_options, loc="upper left")
    grid(True)
    xlabel('degree')
    ylabel(form)
    title('CPU time')

    # Plot speedups
    subplot(122)
    q0, t0 = results[form]["-r quadrature"]
    for (j, test_option) in enumerate(test_options):
        q, t = results[form][test_option]
        t = [t0[k] / t[k] for k in range(len(t))]
        semilogy(q, t, bullets[j])
        hold(True)
    a = list(axis()); a[-1] = 50.0*a[-1]; axis(a);
    legend(test_options, loc="upper left")
    grid(True)
    xlabel('degree')
    title("Speedup vs '-r quadrature'")

show()
