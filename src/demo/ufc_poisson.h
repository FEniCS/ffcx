#include "../ufc/ufc.h"

namespace poisson
{

  class finite_element_0 : public ufc::finite_element
  {
  public:

    inline const char* description() const { return "Lagrange finite element of degree 1 on a triangle"; }

    inline unsigned int space_dimension() const { return 3; }

    inline unsigned int value_rank() const { return 0; }

    inline unsigned int value_dimension(unsigned int i) const { return 0; }

    inline unsigned int num_sub_elements(unsigned int i) const { return 1; }

    inline const finite_element& sub_element(unsigned int i) const { return *this; }

    void evaluate_basis(double* values, const double* x, unsigned int i, const ufc::cell& c) const {}
    
    double evaluate_node(unsigned int n, const ufc::function& f, const ufc::cell& c) const { return 0.0; }

    void tabulate_nodes(unsigned int *nodes, const ufc::mesh& m, const ufc::cell& c) const
    {
      nodes[0] = c.entities[0][0];
      nodes[1] = c.entities[0][1];
      nodes[2] = c.entities[0][2];
    }

    void tabulate_vertex_values(double* vertex_values, const double* nodal_values) const
    {
      vertex_values[0] = nodal_values[0];
      vertex_values[1] = nodal_values[1];
      vertex_values[2] = nodal_values[2];
    }

  };

  class finite_element_1 : public ufc::finite_element
  {
  public:

    inline const char* description() const { return "Lagrange finite element of degree 1 on a triangle"; }

    inline unsigned int space_dimension() const { return 3; }

    inline unsigned int value_rank() const { return 0; }

    inline unsigned int value_dimension(unsigned int i) const { return 0; }

    inline unsigned int num_sub_elements(unsigned int i) const { return 1; }

    inline const finite_element& sub_element(unsigned int i) const { return *this; }

    void evaluate_basis(double* values, const double* x, unsigned int i, const ufc::cell& c) const {}
    
    double evaluate_node(unsigned int n, const ufc::function& f, const ufc::cell& c) const { return 0.0; }

    void tabulate_nodes(unsigned int *nodes, const ufc::mesh& m, const ufc::cell& c) const
    {
      nodes[0] = c.entities[0][0];
      nodes[1] = c.entities[0][1];
      nodes[2] = c.entities[0][2];
    }

    void tabulate_vertex_values(double* vertex_values, const double* nodal_values) const
    {
      vertex_values[0] = nodal_values[0];
      vertex_values[1] = nodal_values[1];
      vertex_values[2] = nodal_values[2];
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

    inline const char* description() const { return "The bilinear form of Poisson's equation"; }

    inline unsigned int rank() const { return 2; }

    inline unsigned int num_coefficients() const { return 0; }

    inline bool interior_contribution() const { return true; }

    inline bool boundary_contribution() const { return false; }

    void tabulate_interior(double* A, const double * const * w, const ufc::cell& c) const
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

    void tabulate_boundary(double* A, const double * const * w, const ufc::cell& c, unsigned int facet) const {}
    
  };

}
