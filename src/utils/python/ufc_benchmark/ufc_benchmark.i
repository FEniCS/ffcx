%module ufc_benchmark

%{
#include "ufc_benchmark.h"
%}

%include "ufc_benchmark.h"

%pythoncode{

def benchmark_forms(forms, geometric_dimension, n=1e6):
    from time import time
    from ufc_benchmark import benchmark
    times = []
    for f in forms:
        tm = 1e300
        for i in range(3):
            t = -time()
            benchmark(f, geometric_dimension, int(n))
            t += time()
            tm = min(t, tm)
        times.append(tm)
    return times

}
