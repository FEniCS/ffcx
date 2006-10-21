#include "ufc.h"

namespace poisson
{

  class node_map_0 : public ufc::node_map
  {
  public:
    
    inline unsigned int space_dimension() const { return 3; }
    
    void tabulate(unsigned int* nodes, const ufc::mesh& m, const ufc::cell& c) const
    {
      nodes[0] = c.entities[0][0];
      nodes[1] = c.entities[0][1];
      nodes[2] = c.entities[0][2];
    }
    
  };

  class node_map_1 : public ufc::node_map
  {
  public:
    
    inline unsigned int space_dimension() const { return 3; }
    
    void tabulate(unsigned int* nodes, const ufc::mesh& m, const ufc::cell& c) const
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
      node_maps = new ufc::node_map * [2];
      node_maps[0] = new node_map_0();
      node_maps[1] = new node_map_1();
    }

    ~element_tensor()
    {
      delete node_maps[0];
      delete node_maps[1];
      delete [] node_maps;
    }

    inline unsigned int rank() const { return 2; }

    inline unsigned int num_coefficients() const { return 0; }

    void tabulate(double* A, const double** w, const ufc::cell& c) const
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
