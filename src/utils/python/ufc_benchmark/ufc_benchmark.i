%module ufc_benchmark

%{
#include "ufc_benchmark.h"
%}

%include "ufc_benchmark.h"

%pythoncode{

def benchmark_forms(forms, geometric_dimension, n=1e8):
    from time import time
    times = []
    for f in forms:
        # take best time of three runs
        ts = 1e999
        for i in range(3):
            t = -time()
            benchmark(f, geometric_dimension, n)
            t += time()
            ts = min(t, ts)
        times.append(ts)
    return times

}
