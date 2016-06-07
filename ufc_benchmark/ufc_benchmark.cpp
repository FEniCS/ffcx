// This is utility code for UFC (Unified Form-assembly Code).
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2006-2015.

#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::vector;

#include <ctime>
#define TMIN 3.0
#define MMIN 1000

#include "ufc_data.h"
#include "ufc_reference_cell.h"
#include "ufc_benchmark.h"

clock_t __tic_time;

void tic()
{
  __tic_time = clock();
}

double toc()
{
  clock_t __toc_time = clock();
  double elapsed_time = ((double) (__toc_time - __tic_time)) / CLOCKS_PER_SEC;
  return elapsed_time;
}

// Adaptive timing: make sure we run for at least TMIN to get reliable results
double time_tabulate_tensor(ufc::cell_integral& integral, double *A, const double * const * w, const ufc::cell & c)
{
  std::size_t M = MMIN;
  while ( true )
  {
    tic();
    for (std::size_t i = 0; i < M; i++)
    {
      integral.tabulate_tensor(A, w, c);
    }
    double t = toc();
    if ( t >= TMIN )
      return t / static_cast<double>(M);
    M *= 10;
    cout << "Elapsed time too short, increasing number of iterations to" << M << endl;
  }

  return 0.0;
}

// Adaptive timing: make sure we run for at least TMIN to get reliable results
double time_tabulate_tensor(ufc::exterior_facet_integral& integral, double *A, const double * const * w, const ufc::cell & c, std::size_t facet)
{
  std::size_t M = MMIN;
  while ( true )
  {
    tic();
    for (std::size_t i = 0; i < M; i++)
    {
      integral.tabulate_tensor(A, w, c, facet);
    }
    double t = toc();
    if ( t >= TMIN )
      return t / static_cast<double>(M);
    M *= 10;
    cout << "Elapsed time too short, increasing number of iterations to" << M << endl;
  }

  return 0.0;
}

// Benchmark all integrals of a form.
vector< vector<double> > benchmark(const ufc::form & form, bool print_tensor)
{
    // construct and allocate some stuff
    ufc::ufc_data data(form);

    // create a reference cell geometry
    ufc::reference_cell c(data.elements[0]->cell_shape());

    // data structures for times
    vector<double> cell_times(form.num_cell_domains());
    vector<double> exterior_facet_times(form.num_exterior_facet_domains());
    vector<double> interior_facet_times(form.num_interior_facet_domains());

    // benchmark all cell integrals
    for(unsigned i = 0; i < form.num_cell_domains(); i++)
    {
        cell_times[i] = time_tabulate_tensor(*data.cell_domains[i], data.A, data.w, c);

        if(print_tensor)
        {
          cout << "Cell element tensor " << i << ":" << endl;
          data.print_tensor();
          cout << endl;
        }
    }

    // benchmark all exterior facet integrals
    for(unsigned i = 0; i < form.num_exterior_facet_domains(); i++)
    {
        std::size_t facet = 0; // TODO: would it be interesting to time all facets?
        exterior_facet_times[i] = time_tabulate_tensor(*data.exterior_facet_domains[i], data.A, data.w, c, facet);

        if(print_tensor)
        {
          cout << "Exterior facet element tensor " << i << ":" << endl;
          data.print_tensor();
          cout << endl;
        }
    }

    // benchmark all interior facet integrals
    /* // TODO: If somebody needs this, please implement it! Need two cells, and larger A.
    for(unsigned i = 0; i < form.num_interior_facet_domains(); i++)
    {
        std::size_t facet = 0; // TODO: would it be interesting to time all facets?
        interior_facet_times[i] = time_tabulate_tensor(*data.interior_facet_domains[i], data.A, data.w, c, facet);

        if(print_tensor)
        {
          cout << "Interior facet element tensor " << i << ":" << endl;
          data.print_tensor();
          cout << endl;
        }
    }
    */

    vector< vector<double> > result(3);
    result[0] = cell_times;
    result[1] = exterior_facet_times;
    result[2] = interior_facet_times;

    return result;
}


vector< vector<double> > tabulate_cell_tensor(const ufc::form & form, vector< vector<double> > w, int domain)
{
  ufc::ufc_data data(form);

  // copy w to the appropriate array
  if(data.num_coefficients != w.size())
      throw std::runtime_error("Wrong number of coefficients");
  for(unsigned i=0; i<data.num_coefficients; i++)
  {
    if(data.dimensions[data.rank+i] != w[i].size())
        throw std::runtime_error("Wrong coefficient dimension.");
    for(unsigned j=0; j<data.dimensions[data.rank+i]; j++)
    {
      data.w[i][j] = w[i][j];
    }
  }

  // create a reference cell geometry
  ufc::reference_cell c(data.elements[0]->cell_shape());

  // tabulate the tensor
  data.cell_integrals[domain]->tabulate_tensor(data.A, data.w, c);

  // copy element tensor to stl-structure for easy returning to python (should perhaps rather use numpy and some typemaps, but I'm lazy)
  vector< vector<double> > A;
  if(data.rank == 2)
  {
    A.resize(data.dimensions[0]);
    for(unsigned i=0; i<data.dimensions[0]; i++)
    {
      A[i].resize(data.dimensions[1]);
      for(unsigned j=0; j<data.dimensions[1]; j++)
      {
        A[i][j] = data.A[i*data.dimensions[1] + j];
      }
    }
  }
  else if(data.rank == 1)
  {
    A.resize(data.dimensions[0]);
    for(unsigned i=0; i<data.dimensions[0]; i++)
    {
      A[i].resize(1);
      A[i][0] = data.A[i];
    }
  }
  else if(data.rank == 0)
  {
    A.resize(1);
    A[0].resize(1);
    A[0][0] = data.A[0];
  }
  else
  {
    throw std::runtime_error("rank != 0,1,2 not implemented");
  }

  return A;
}

std::vector< std::vector<double> > tabulate_cell_integral(const std::shared_ptr<ufc::form> form, std::vector< std::vector<double> > w, ufc::cell cell, int domain)
{
  ufc::ufc_data data(*form);

  // copy w to the appropriate array
  if(data.num_coefficients != w.size())
      throw std::runtime_error("Wrong number of coefficients");
  for(unsigned i=0; i<data.num_coefficients; i++)
  {
    if(data.dimensions[data.rank+i] != w[i].size())
        throw std::runtime_error("Wrong coefficient dimension.");
    for(unsigned j=0; j<data.dimensions[data.rank+i]; j++)
    {
      data.w[i][j] = w[i][j];
    }
  }

  // tabulate the tensor
  data.cell_integrals[domain]->tabulate_tensor(data.A, data.w, cell);

  // copy element tensor to stl-structure for easy returning to python (should perhaps rather use numpy and some typemaps, but I'm lazy)
  vector< vector<double> > A;
  if(data.rank == 2)
  {
    A.resize(data.dimensions[0]);
    for(unsigned i=0; i<data.dimensions[0]; i++)
    {
      A[i].resize(data.dimensions[1]);
      for(unsigned j=0; j<data.dimensions[1]; j++)
      {
        A[i][j] = data.A[i*data.dimensions[1] + j];
      }
    }
  }
  else if(data.rank == 1)
  {
    A.resize(data.dimensions[0]);
    for(unsigned i=0; i<data.dimensions[0]; i++)
    {
      A[i].resize(1);
      A[i][0] = data.A[i];
    }
  }
  else if(data.rank == 0)
  {
    A.resize(1);
    A[0].resize(1);
    A[0][0] = data.A[0];
  }
  else
  {
    throw std::runtime_error("rank != 0,1,2 not implemented");
  }

  return A;
}

std::vector< std::vector<double> > tabulate_exterior_facet_integral(const std::shared_ptr<ufc::form> form, std::vector< std::vector<double> > w, ufc::cell& cell, int facet, int domain)
{
  ufc::ufc_data data(*form);

  // copy w to the appropriate array
  if(data.num_coefficients != w.size())
      throw std::runtime_error("Wrong number of coefficients");
  for(unsigned i=0; i<data.num_coefficients; i++)
  {
    if(data.dimensions[data.rank+i] != w[i].size())
        throw std::runtime_error("Wrong coefficient dimension.");
    for(unsigned j=0; j<data.dimensions[data.rank+i]; j++)
    {
      data.w[i][j] = w[i][j];
    }
  }

  // tabulate the tensor
  data.exterior_facet_integrals[domain]->tabulate_tensor(data.A, data.w, cell, facet);

  // copy element tensor to stl-structure for easy returning to python (should perhaps rather use numpy and some typemaps, but I'm lazy)
  vector< vector<double> > A;
  if(data.rank == 2)
  {
    A.resize(data.dimensions[0]);
    for(unsigned i=0; i<data.dimensions[0]; i++)
    {
      A[i].resize(data.dimensions[1]);
      for(unsigned j=0; j<data.dimensions[1]; j++)
      {
        A[i][j] = data.A[i*data.dimensions[1] + j];
      }
    }
  }
  else if(data.rank == 1)
  {
    A.resize(data.dimensions[0]);
    for(unsigned i=0; i<data.dimensions[0]; i++)
    {
      A[i].resize(1);
      A[i][0] = data.A[i];
    }
  }
  else if(data.rank == 0)
  {
    A.resize(1);
    A[0].resize(1);
    A[0][0] = data.A[0];
  }
  else
  {
    throw std::runtime_error("rank != 0,1,2 not implemented");
  }

  return A;
}

std::vector< std::vector<double> > tabulate_interior_facet_integral(const std::shared_ptr<ufc::form> form, std::vector< std::vector<double> > macro_w,\
                                                                    ufc::cell& cell0, ufc::cell& cell1, int facet_0, int facet_1, int domain)
{
  ufc::ufc_data data(*form);

  // copy w to the appropriate array
  if(data.num_coefficients != macro_w.size())
      throw std::runtime_error("Wrong number of coefficients");
  for(unsigned i=0; i<data.num_coefficients; i++)
  {
    if(2*data.dimensions[data.rank+i] != macro_w[i].size())
        throw std::runtime_error("Wrong coefficient dimension.");
    for(unsigned j=0; j<2*data.dimensions[data.rank+i]; j++)
    {
      data.macro_w[i][j] = macro_w[i][j];
    }
  }

  // tabulate the tensor
  data.interior_facet_integrals[domain]->tabulate_tensor(data.macro_A, data.macro_w, cell0, cell1, facet_0, facet_1);

  // copy element tensor to stl-structure for easy returning to python (should perhaps rather use numpy and some typemaps, but I'm lazy)
  vector< vector<double> > A;
  if(data.rank == 2)
  {
    A.resize(2*data.dimensions[0]);
    for(unsigned i=0; i<2*data.dimensions[0]; i++)
    {
      A[i].resize(2*data.dimensions[1]);
      for(unsigned j=0; j<2*data.dimensions[1]; j++)
      {
        A[i][j] = data.macro_A[i*2*data.dimensions[1] + j];
      }
    }
  }
  else if(data.rank == 1)
  {
    A.resize(2*data.dimensions[0]);
    for(unsigned i=0; i<2*data.dimensions[0]; i++)
    {
      A[i].resize(1);
      A[i][0] = data.macro_A[i];
    }
  }
  else if(data.rank == 0)
  {
    A.resize(1);
    A[0].resize(1);
    A[0][0] = data.macro_A[0];
  }
  else
  {
    throw std::runtime_error("rank != 0,1,2 not implemented");
  }

  return A;
}
