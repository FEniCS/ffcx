// This is utility code for UFC (Unified Form-assembly Code).
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2006-2015.

#ifndef __UFC_BENCHMARK_H__
#define __UFC_BENCHMARK_H__

#include "ufc.h"
#include <memory>
#include <vector>

/* Benchmark time to run tabulate_tensor for all integrals in a form. *
 * Uses a reference cell and one-cell mesh, and sets all w_ij = 1.0.  */
std::vector< std::vector<double> > benchmark(const ufc::form & form,
                                             bool print_tensors);

/* Compute one element tensor on the reference cell with the given coefficients. */
std::vector< std::vector<double> > tabulate_cell_tensor(const ufc::form & form,
                                                        std::vector< std::vector<double> > w,
                                                        int domain);

/* Compute one cell integral. */
std::vector< std::vector<double> > tabulate_cell_integral(const std::shared_ptr<ufc::form> form,
                                                          std::vector< std::vector<double> > w,
                                                          ufc::cell cell,
                                                          int domain);

/* Compute one exterior facet integral. */
std::vector< std::vector<double> > tabulate_exterior_facet_integral(const std::shared_ptr<ufc::form> form,
                                                                    std::vector< std::vector<double> > w,
                                                                    ufc::cell& cell,
                                                                    int facet,
                                                                    int domain);

/* Compute one interior facet integral. */
std::vector< std::vector<double> > tabulate_interior_facet_integral(const std::shared_ptr<ufc::form> form,
                                                                    std::vector< std::vector<double> > macro_w,
                                                                    ufc::cell& cell0,
                                                                    ufc::cell& cell1,
                                                                    int facet_0,
                                                                    int facet_1,
                                                                    int domain);

#endif
