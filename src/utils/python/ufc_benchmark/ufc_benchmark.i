%module ufc_benchmark

%{
#include "ufc_benchmark.h"
%}

%include "ufc_benchmark.h"

%pythoncode{

def benchmark_forms(forms, geometric_dimension, n=1e8):
    from timeit import Timer
    for f in forms:
        t = Timer( "benchmark(f, geometric_dimension, n)" )
        times += [ min( t.repeat(3, 1) ) ]
    return times

}
