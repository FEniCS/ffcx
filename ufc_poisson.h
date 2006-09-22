#include "ufc.h"

namespace poisson
{

  class finite_element_0 : public ufc::finite_element
  {
  public:
    
    inline unsigned int space_dimension() const { return 3; }
    
    inline unsigned int shape_dimension() const { return 2; }

    void compute_node_map() const {}
    
  };

  class finite_element_1 : public ufc::finite_element
  {
  public:
    
    inline unsigned int space_dimension() const { return 3; }
    
    inline unsigned int shape_dimension() const { return 2; }
    
    void compute_node_map() const {}

  };
  
  class multilinear_form : public ufc::multilinear_form
  {
  public:
    
    multilinear_form() : ufc::multilinear_form()
    {
      finite_elements = new ufc::finite_element * [2];
      finite_elements[0] = new finite_element_0();
      finite_elements[1] = new finite_element_1();
    }

    ~multilinear_form()
    {
      delete finite_elements[0];
      delete finite_elements[1];
      delete [] finite_elements;
    }

    inline unsigned int num_arguments() const { return 2; }

    inline unsigned int num_coefficients() const { return 0; }
    
    void compute_element_tensor_interior(double* A) const
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
    
    void compute_element_tensor_boundary(double* A) const
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

  };

}
