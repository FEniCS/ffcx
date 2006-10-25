#include "ufc.h"

namespace poisson
{

  class finite_element_0 : public ufc::finite_element
  {
  public:

    inline unsigned int space_dimension() const { return 3; }

    double evaluate_basis(unsigned int n, double* x, const ufc::cell& c) const { return 0.0; }
    
    double evaluate_node(unsigned int n, const ufc::function& f, const ufc::cell& c) const { return 0.0; }

    void tabulate_nodes(unsigned int *nodes, const ufc::mesh& m, const ufc::cell& c) const
    {
      nodes[0] = c.entities[0][0];
      nodes[1] = c.entities[0][1];
      nodes[2] = c.entities[0][2];
    }

  };

  class finite_element_1 : public ufc::finite_element
  {
  public:

    inline unsigned int space_dimension() const { return 3; }

    double evaluate_basis(unsigned int n, double* x, const ufc::cell& c) const { return 0.0; }
    
    double evaluate_node(unsigned int n, const ufc::function& f, const ufc::cell& c) const { return 0.0; }

    void tabulate_nodes(unsigned int *nodes, const ufc::mesh& m, const ufc::cell& c) const
    {
      nodes[0] = c.entities[0][0];
      nodes[1] = c.entities[0][1];
      nodes[2] = c.entities[0][2];
    }

  };
  
  class element_tensor : public ufc::element_tensor
  {
  public:
    
    element_tensor() : ufc::element_tensor()
    {
      finite_elements = new ufc::finite_element * [2];
      finite_elements[0] = new finite_element_0();
      finite_elements[1] = new finite_element_1();
    }

    ~element_tensor()
    {
      delete finite_elements[0];
      delete finite_elements[1];
      delete [] finite_elements;
    }

    inline unsigned int rank() const { return 2; }

    inline unsigned int num_coefficients() const { return 0; }

    inline bool interior_contribution() const { return true; }

    inline bool boundary_contribution() const { return false; }

    void tabulate_interior(double* A, const double** w, const ufc::cell& c) const
    {
      A[0] = 0.0;
      A[1] = 0.0;
      A[2] = 0.0;
      A[3] = 0.0;
      A[4] = 0.0;
      A[5] = 0.0;
      A[6] = 0.0;
      A[7] = 0.0;
      A[8] = 0.0;
    }

    void tabulate_boundary(double* A, const double** w, const ufc::cell& c, unsigned int facet) const {}
    
  };

}
