#include "ufc_benchmark.h"

typedef unsigned int uint; 

int benchmark(const ufc::form & form, uint geometric_dimension, uint n)
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
        	w[i][j] = 1.0;
        }
    }

    // TODO: assuming one cell integral only, allow benchmarking of any single integral
    ufc::cell_integral *itg = form.create_cell_integral(0);

    // assuming top. and geom. dims are equal
    ufc::mesh m;
    m.topological_dimension = geometric_dimension;
    m.geometric_dimension   = geometric_dimension;

    ufc::cell c;
    c.topological_dimension = m.topological_dimension;
    c.geometric_dimension   = m.geometric_dimension;

    ufc::finite_element * fe = form.create_finite_element(0);
    c.cell_shape = fe->cell_shape();
    delete fe;

    // TODO: perhaps create utility cell classes containing the
    //       standard entity numberings and reference geometries?
    uint num_vertices = 0;
    switch(c.cell_shape)
    {
    case ufc::interval:
    	num_vertices = 2;
    	break;
    case ufc::triangle:
    	num_vertices = 3;
    	break;
    case ufc::tetrahedron:
    	num_vertices = 4;
    	break;
    case ufc::quadrilateral:
    	num_vertices = 4;
    	break;
    case ufc::hexahedron:
    	num_vertices = 8;
    	break;
    }
    c.coordinates = new double*[num_vertices];
    for(uint i=0; i<num_vertices; i++)
    {
    	c.coordinates[i] = new double[geometric_dimension];
    }
    // FIXME: fill coordinates with something sensible
    
    
    // benchmark!
    for(unsigned int i=0; i<n; i++)
    {
        itg->tabulate_tensor(A, w, c);
    }

    // clean up
    delete itg;
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

    for(uint i=0; i<num_vertices; i++)
    {
    	delete [] c.coordinates[i];
    }
    delete [] c.coordinates;
}

/*
#define TMIN 1.0
#define MMIN 10

#include <ctime>

clock_t __tic_time;

void tic()
{
  __tic_time = clock();
}

double toc()
{
  clock_t __toc_time = clock();
  double elapsed_time = ((real) (__toc_time - __tic_time)) / CLOCKS_PER_SEC;
  return elapsed_time;
}

// Adaptive timing: make sure we run for at least TMIN to get reliable results
double time_tabulate_tensor(ufc::cell_integral& cell_integral)
{
  unsigned int M = MMIN;
  while ( true )
  {
    tic();
    for (unsigned int i = 0; i < M; i++)
    {
      cell_integral.tabulate_tensor(A, w, c); // FIXME: define A, w, c
    }
    double t = toc();
    if ( t >= TMIN )
      return t / static_cast<double>(M);
    M *= 10;
    cout << "Elapsed time too short, increasing number of iterations to" << M << endl;
  }

  return 0.0;
}

*/
