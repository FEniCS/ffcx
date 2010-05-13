"This script plots the results found in bench.log."

__author__ = "Anders Logg"
__date__ = "2010-05-13"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from pylab import *

# Read logfile
results = {}
for line in open("bench.log", "r").read().split("\n"):
    if not "," in line: continue
    test_case, test_option, timing = [w.strip() for w in line.split(",")]
    form, degree = test_case.split("_")
    if not form in results:
        results[form] = {}
    if not test_option in results[form]:
        results[form][test_option] = ([], [])
    results[form][test_option][0].append(int(degree))
    results[form][test_option][1].append(timing)

# Plot results
forms = sorted([form for form in results])
test_options = ["-r tensor", "-r tensor -O", "-r quadrature", "-r quadrature -O"]
bullets = ["x-", "o-", "*-", "s-"]
for (i, form) in enumerate(forms):
    figure(i)
    for (j, test_option) in enumerate(test_options):
        q, t = results[form][test_option]
        semilogy(q, t, bullets[j])
        hold(True)
    legend(test_options, loc="upper left")
    grid(True)
    title(form)
    xlabel('degree')
    ylabel('CPU time')

show()
