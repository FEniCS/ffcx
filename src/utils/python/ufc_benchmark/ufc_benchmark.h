// This is utility code for UFC (Unified Form-assembly Code) v. 1.0.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenics.org/) 2006-2007.

#ifndef __UFC_BENCHMARK_H__
#define __UFC_BENCHMARK_H__

#include <ufc.h>
#include <vector>

/*
 * Benchmark time to run tabulate_tensor for all integrals in a form.
 * Uses a reference cell and one-cell mesh, and sets all w_ij = 1.0.
 */
std::vector< std::vector<double> > benchmark(const ufc::form & form, bool print_tensors);

/*
 * Compute one element tensor on the reference cell with the given coefficients.
 */
std::vector< std::vector<double> > tabulate_cell_tensor(const ufc::form & form, std::vector< std::vector<double> > w, int domain);

#endif

