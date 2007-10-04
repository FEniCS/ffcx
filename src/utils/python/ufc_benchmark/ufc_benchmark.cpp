

#define TMIN 1.0
#define MMIN 10

#include <iostream>
using std::cout;
using std::endl;

#include <ctime>

#include "ufc_benchmark.h"
#include "ufc_reference_cell.h"

typedef unsigned int uint; 


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
  unsigned int M = MMIN;
  while ( true )
  {
    tic();
    for (unsigned int i = 0; i < M; i++)
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
double time_tabulate_tensor(ufc::exterior_facet_integral& integral, double *A, const double * const * w, const ufc::cell & c, unsigned int facet)
{
  unsigned int M = MMIN;
  while ( true )
  {
    tic();
    for (unsigned int i = 0; i < M; i++)
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
std::vector< std::list<double> > benchmark(const ufc::form & form)
{
    uint rank  = form.rank();
    uint num_w = form.num_coefficients();

    // construct dofmaps to get dimensions
    ufc::dof_map ** dofmaps = new ufc::dof_map*[rank+num_w];
    for(uint i=0; i<rank+num_w; i++)
    {
        dofmaps[i] = form.create_dof_map(i);
    }

    // allocate A
    uint A_size = 1;
    for(uint i=0; i<rank; i++)
    {
        A_size *= dofmaps[i]->local_dimension();
    }
    double* A = new double[A_size];

    // allocate local coefficient data
    double ** w = new double*[num_w];
    for(uint i=0; i<num_w; i++)
    {
        uint dim = dofmaps[i+rank]->local_dimension();
        w[i] = new double[dim];
        for(uint j=0; j<dim; j++)
        {
        	// WARNING: assuming this is valid input data... Could potentially lead to division by zero and such.
        	w[i][j] = 0.0;
        }
    }

    // create a reference cell geometry
    ufc::finite_element * fe = form.create_finite_element(0);
    ufc::reference_cell c(fe->cell_shape());
    delete fe;

    // data structures for times
    std::list<double> cell_times;
    std::list<double> exterior_facet_times;
    std::list<double> interior_facet_times;

    // benchmark all cell integrals
    for(uint i = 0; i < form.num_cell_integrals(); i++)
    {
        ufc::cell_integral *itg = form.create_cell_integral(i);
        double t = time_tabulate_tensor(*itg, A, w, c);
        delete itg;
        cell_times.push_back(t);
    }

    // benchmark all exterior facet integrals
    for(uint i = 0; i < form.num_exterior_facet_integrals(); i++)
    {
        ufc::exterior_facet_integral *itg = form.create_exterior_facet_integral(i);
        unsigned int facet = 0; // TODO: would it be interesting to time all facets?
        double t = time_tabulate_tensor(*itg, A, w, c, facet);
        delete itg;
        exterior_facet_times.push_back(t);
    }

    // benchmark all interior facet integrals
    /* // TODO: If somebody needs this, please implement it! Need two cells, and larger A.
    for(uint i = 0; i < form.num_interior_facet_integrals(); i++)
    {
        ufc::interior_facet_integral *itg = form.create_interior_facet_integral(i);
        double t = time_tabulate_tensor(*itg, A, w, c);
        delete itg;
        interior_facet_times.push_back(t);
    }
    */

    std::vector< std::list<double> > result(3);
    result[0] = cell_times;
    result[1] = exterior_facet_times;
    result[2] = interior_facet_times;


    // clean up
    delete [] A;

    for(uint i=0; i<num_w; i++)
    {
        delete [] w[i];
    }
    delete [] w;

    for(uint i=0; i<rank+num_w; i++)
    {
        delete dofmaps[i];
    }
    delete [] dofmaps;


    return result;
}
