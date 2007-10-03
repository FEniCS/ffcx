%module ufc_benchmark

%{
#include "ufc_benchmark.h"
%}

%include "ufc_benchmark.h"

%include stl.i
%include std_list.i
%include std_vector.i


%template(list_double) std::list<double>;
%template(vector_list_double) std::vector< std::list<double> >;

%typedef std::list<double> list_double;
%typedef std::vector< std::list<double> > vector_list_double;


%pythoncode{

def benchmark_forms(forms):
    import gc
    gc.collect()
    
    times = []
    for f in forms:
        res = benchmark(f)
        times.append(res)
    return times

}
