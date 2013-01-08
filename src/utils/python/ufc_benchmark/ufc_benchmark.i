%module ufc_benchmark

// ------------------------ STL stuff

%{
#include <vector>
%}

%include stl.i
%include std_vector.i
%include std_carray.i

%template(vector_double)     std::vector<double>;
%typedef std::vector<double> vector_double;

%template(vector_vector_double)             std::vector< std::vector<double> >;
%typedef std::vector< std::vector<double> > vector_vector_double;

%template(vector_std_t) std::vector< std::size_t >;
%typedef std::vector< std::size_t > vector_size_t;

%template(vector_vector_size_t) std::vector< std::vector< std::size_t > >;
%typedef std::vector< std::vector< std::size_t > > vector_vector_size_t;

// ------------------------ UFC stuff

%import ufc.i

%{
#include "ufc.h"
#include "ufc_benchmark.h"
#include "ufc_reference_cell.h"
%}

%include "ufc.h"
%include "ufc_benchmark.h"
%include "ufc_reference_cell.h"

// ----------------------- Reference to shared pointer utility

%{
class NoDeleter { public: void operator()(ufc::form *) {} };
boost::shared_ptr<ufc::form> form_ptr(ufc::form * form) { return boost::shared_ptr<ufc::form>(form, NoDeleter()); }
%}
class NoDeleter { public: void operator()(ufc::form *) {} };
boost::shared_ptr<ufc::form> form_ptr(ufc::form * form) { return boost::shared_ptr<ufc::form>(form, NoDeleter()); }

// ----------------------- Python wrapper for benchmark

%pythoncode{

def benchmark_forms(forms, print_tensors):
    import gc
    gc.collect()

    times = []
    for f in forms:
        res = benchmark(f, print_tensors)
        times.append(tuple(res))
    return times

}
