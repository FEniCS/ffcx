// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.3-dev.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __PROJECTION_H
#define __PROJECTION_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/FiniteElementSpec.h>
#include <dolfin/BilinearForm.h>
#include <dolfin/LinearForm.h>
#include <dolfin/Functional.h>
#include <dolfin/FEM.h>

namespace dolfin { namespace Projection {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class LinearForm : public dolfin::LinearForm
{
public:

  class TestElement;

  class FunctionElement_0;

  LinearForm(Function& w0);
  

  void eval(real block[], const AffineMap& map) const;

  void eval(real block[], const AffineMap& map, unsigned int facet) const;

};

class LinearForm::TestElement : public dolfin::FiniteElement
{
public:

  TestElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~TestElement()
  {
    if ( tensordims ) delete [] tensordims;
    if ( subelements )
    {
      for (unsigned int i = 0; i < elementdim(); i++)
        delete subelements[i];
      delete [] subelements;
    }
  }

  inline unsigned int spacedim() const
  {
    return 3;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    dolfin_error("Element is scalar.");
    return 0;
  }

  inline unsigned int elementdim() const
  {
    return 1;
  }

  inline unsigned int rank() const
  {
    return 0;
  }

  void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
  {
    nodes[0] = cell.vertexID(0);
    nodes[1] = cell.vertexID(1);
    nodes[2] = cell.vertexID(2);
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
  }

  const FiniteElement& operator[] (unsigned int i) const
  {
    return *this;
  }

  FiniteElement& operator[] (unsigned int i)
  {
    return *this;
  }

  FiniteElementSpec spec() const
  {
    FiniteElementSpec s("Lagrange", "triangle", 1);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class LinearForm::FunctionElement_0 : public dolfin::FiniteElement
{
public:

  FunctionElement_0() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_0()
  {
    if ( tensordims ) delete [] tensordims;
    if ( subelements )
    {
      for (unsigned int i = 0; i < elementdim(); i++)
        delete subelements[i];
      delete [] subelements;
    }
  }

  inline unsigned int spacedim() const
  {
    return 3;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    dolfin_error("Element is scalar.");
    return 0;
  }

  inline unsigned int elementdim() const
  {
    return 1;
  }

  inline unsigned int rank() const
  {
    return 0;
  }

  void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
  {
    nodes[0] = cell.vertexID(0);
    nodes[1] = cell.vertexID(1);
    nodes[2] = cell.vertexID(2);
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
  }

  const FiniteElement& operator[] (unsigned int i) const
  {
    return *this;
  }

  FiniteElement& operator[] (unsigned int i)
  {
    return *this;
  }

  FiniteElementSpec spec() const
  {
    FiniteElementSpec s("Lagrange", "triangle", 1);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

LinearForm::LinearForm(Function& w0) : dolfin::LinearForm(1)
{
  // Create finite element for test space
  _test = new TestElement();

  // Add functions
  initFunction(0, w0, new FunctionElement_0());
}

void LinearForm::eval(real block[], const AffineMap& map) const
{
  // Compute coefficients
  const real c0_0 = 3.333333333333333e-01*c[0][0] + 3.333333333333334e-01*c[0][1] + 3.333333333333333e-01*c[0][2];
  const real c1_0 = c[0][0];
  const real c1_1 = c[0][1];
  const real c1_2 = c[0][2];
  const real c2_0 = c[0][0];
  const real c2_1 = c[0][1];
  const real c2_2 = c[0][2];
  const real c2_3 = 5.000000000000000e-01*c[0][1] + 4.999999999999998e-01*c[0][2];
  const real c2_4 = 4.999999999999999e-01*c[0][0] + 4.999999999999997e-01*c[0][2];
  const real c2_5 = 5.000000000000000e-01*c[0][0] + 5.000000000000001e-01*c[0][1];

  // Compute geometry tensors
  const real G0_0 = map.det*c0_0;
  const real G1_0 = map.det*c1_0;
  const real G1_1 = map.det*c1_1;
  const real G1_2 = map.det*c1_2;
  const real G2_0 = map.det*c2_0;
  const real G2_1 = map.det*c2_1;
  const real G2_2 = map.det*c2_2;
  const real G2_3 = map.det*c2_3;
  const real G2_4 = map.det*c2_4;
  const real G2_5 = map.det*c2_5;

  // Compute element tensor
  block[0] = 1.666666666666665e-01*G0_0 + 8.333333333333318e-02*G1_0 + 4.166666666666659e-02*G1_1 + 4.166666666666658e-02*G1_2 + 1.666666666666664e-02*G2_0 - 8.333333333333316e-03*G2_1 - 8.333333333333316e-03*G2_2 + 3.333333333333327e-02*G2_3 + 6.666666666666654e-02*G2_4 + 6.666666666666654e-02*G2_5;
  block[1] = 1.666666666666665e-01*G0_0 + 4.166666666666659e-02*G1_0 + 8.333333333333318e-02*G1_1 + 4.166666666666659e-02*G1_2 - 8.333333333333318e-03*G2_0 + 1.666666666666664e-02*G2_1 - 8.333333333333312e-03*G2_2 + 6.666666666666654e-02*G2_3 + 3.333333333333327e-02*G2_4 + 6.666666666666654e-02*G2_5;
  block[2] = 1.666666666666665e-01*G0_0 + 4.166666666666658e-02*G1_0 + 4.166666666666659e-02*G1_1 + 8.333333333333316e-02*G1_2 - 8.333333333333319e-03*G2_0 - 8.333333333333316e-03*G2_1 + 1.666666666666663e-02*G2_2 + 6.666666666666654e-02*G2_3 + 6.666666666666654e-02*G2_4 + 3.333333333333326e-02*G2_5;
}

// No contribution from the boundary
void LinearForm::eval(real block[], const AffineMap& map, unsigned int facet) const {}   
} }

#endif
