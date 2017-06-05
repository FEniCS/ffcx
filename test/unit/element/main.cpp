#include <gtest/gtest.h>

#include <ufc.h>
#include <ufc_geometry.h>
#include <boost/multi_array.hpp>
#include "lagrange_elements.h"
//#include "mock_cells.h"
//#include "debugging.h"

#include "../../uflacs/crosslanguage/cppsupport/mock_cells.h"

//#include "test_ufc_integral_types.h"


namespace
{
  // Tests computes derivatives on a reference cell (using reference
  // and physical version)
  void test_reference_derivatives(const ufc::finite_element& e,
                                  const boost::multi_array<double, 2>& X,
                                  const boost::multi_array<double, 4>& v_ref,
                                  const mock_cell& cell)
  {
    // Number of points
    const std::size_t num_points = X.shape()[0];

    const std::size_t gdim = e.geometric_dimension();
    const std::size_t space_dim = e.space_dimension();
    const std::size_t ref_value_size = e.reference_value_size();

    ASSERT_EQ(v_ref.shape()[0], num_points);
    ASSERT_EQ(v_ref.shape()[1], space_dim);
    ASSERT_EQ(v_ref.shape()[2], gdim);
    ASSERT_EQ(v_ref.shape()[3], ref_value_size);

    // Compute derivatives on reference element
    boost::multi_array<double, 4> v = v_ref;
    std::fill(v.data(), v.data() + v.num_elements(), -10.0);
    e.evaluate_reference_basis_derivatives(v.data(), 1, num_points, X.data());
    for (std::size_t point = 0; point < num_points; ++point)
    {
      for (std::size_t basis = 0; basis < space_dim; ++basis)
      {
        for (std::size_t comp = 0; comp < gdim; ++comp)
        {
          for (std::size_t k = 0; k < ref_value_size; ++k)
            EXPECT_FLOAT_EQ(v_ref[point][basis][comp][k], v[point][basis][comp][k]);
        }
      }
    }

   // Compute derivatives on real element using mock element which
   // corresponds to reference cell
    boost::multi_array<double, 3> w(boost::extents[space_dim][gdim][ref_value_size]);
    for (std::size_t p = 0; p < num_points; ++p)
    {
      std::vector<double> _x(X[p].begin(), X[p].end());
      e.evaluate_basis_derivatives_all(1, w.data(), _x.data(),
                                       cell.coordinate_dofs, 0);
      for (std::size_t basis = 0; basis < space_dim; ++basis)
      {
        for (std::size_t comp = 0; comp < gdim; ++comp)
        {
          for (std::size_t k = 0; k < ref_value_size; ++k)
            EXPECT_NEAR(v_ref[p][basis][comp][k], w[basis][comp][k], 1.0e-13);
        }
      }
    }
  }
}


TEST(ScalarLagrangeInterval, cell_shape)
{
  interval_1_finite_element_0 e1_int;
  EXPECT_EQ(e1_int.cell_shape(), ufc::shape::interval);
}

TEST(ScalarLagrangeIntervalP1, basis_derivatives)
{
  // Create element
  interval_1_finite_element_0 e;

  // dims
  const std::size_t gdim = e.geometric_dimension();
  const std::size_t space_dim = e.space_dimension();
  const std::size_t ref_value_size = e.reference_value_size();

  // Points on reference element to test
  boost::multi_array<double, 2> X(boost::extents[3][1]);
  X[0][0] = -0.2; X[1][0] = 0.5; X[2][0] = -0.1;

  // Reference solution
  boost::multi_array<double, 4> v_ref(boost::extents[3][space_dim][gdim][ref_value_size]);
  v_ref[0][0][0][0] = -1.0;
  v_ref[0][1][0][0] = 1.0;
  v_ref[1][0][0][0] = -1.0;
  v_ref[1][1][0][0] = 1.0;
  v_ref[2][0][0][0] = -1.0;
  v_ref[2][1][0][0] = 1.0;

  mock_cell cell;
  cell.fill_reference_interval(1);
  test_reference_derivatives(e, X, v_ref, cell);
}

TEST(ScalarLagrangeTriangleP1, basis_derivatives)
{
  // Create element
  triangle_1_finite_element_0 e;

  // dims
  const std::size_t gdim = e.geometric_dimension();
  const std::size_t space_dim = e.space_dimension();
  const std::size_t ref_value_size = e.reference_value_size();

  // Points on reference element to test
  boost::multi_array<double, 2> X(boost::extents[3][2]);
  X[0][0] = 0.0;  X[0][1] = 0.0;
  X[1][0] = 0.25; X[1][1] = 0.25;
  X[2][0] = 0.9;  X[2][0] = 0.01;

  // Reference solution
  boost::multi_array<double, 4> v_ref(boost::extents[3][space_dim][gdim][ref_value_size]);
  for (std::size_t p = 0; p < X.shape()[0]; ++p)
  {
    v_ref[p][0][0][0] = -1.0;
    v_ref[p][0][1][0] = -1.0;

    v_ref[p][1][0][0] = 1.0;
    v_ref[p][1][1][0] = 0.0;

    v_ref[p][2][0][0] = 0.0;
    v_ref[p][2][1][0] = 1.0;
  }

  mock_cell cell;
  cell.fill_reference_triangle(2);
  test_reference_derivatives(e, X, v_ref, cell);
}


/*
TEST(ScalarLagrangeInterval, cell_shape)
{
  interval_1_finite_element_0 e1_int;
  EXPECT_EQ(e1_int.cell_shape(), ufc::shape::interval);

  triangle_1_finite_element_0 e1_tri;
  EXPECT_EQ(e1_tri.cell_shape(), ufc::shape::triangle);

  tetrahedron_1_finite_element_0 e1_tet;
  EXPECT_EQ(e1_tet.cell_shape(), ufc::shape::tetrahedron);
}


TEST(ScalarLagrangeInterval, cell_shape)
{
  interval_1_finite_element_0 e1_int;
  EXPECT_EQ(e1_int.cell_shape(), ufc::shape::interval);

  triangle_1_finite_element_0 e1_tri;
  EXPECT_EQ(e1_tri.cell_shape(), ufc::shape::triangle);

  tetrahedron_1_finite_element_0 e1_tet;
  EXPECT_EQ(e1_tet.cell_shape(), ufc::shape::tetrahedron);
}



TEST(P1_Lagrange, cell_shape)
{
  interval_1_finite_element_0 e1_int;
  EXPECT_EQ(e1_int.cell_shape(), ufc::shape::interval);

  triangle_1_finite_element_0 e1_tri;
  EXPECT_EQ(e1_tri.cell_shape(), ufc::shape::triangle);

  tetrahedron_1_finite_element_0 e1_tet;
  EXPECT_EQ(e1_tet.cell_shape(), ufc::shape::tetrahedron);
}

TEST(P1_Lagrange, evaluate_reference_basis_derivatives)
{
  {
    interval_1_finite_element_0 e;
    std::vector<double> X = {0.2, 0.5, 1.0};
    std::vector<double> v(X.size()*2);
    e.evaluate_reference_basis_derivatives(v.data(), 1, X.size(), X.data());
    EXPECT_FLOAT_EQ(v[0], -1.0);
    EXPECT_FLOAT_EQ(v[1], 1.0);

    EXPECT_FLOAT_EQ(v[2], -1.0);
    EXPECT_FLOAT_EQ(v[3], 1.0);

    EXPECT_FLOAT_EQ(v[4], -1.0);
    EXPECT_FLOAT_EQ(v[5], 1.0);
  }

  {
    triangle_1_finite_element_0 e;
    std::vector<double> X = {0.2, 0.5};
    std::vector<double> v((X.size()/2)*3*2);
    e.evaluate_reference_basis_derivatives(v.data(), 1, X.size()/2, X.data());
    EXPECT_FLOAT_EQ(v[0], -1.0);
    EXPECT_FLOAT_EQ(v[1], -1.0);

    EXPECT_FLOAT_EQ(v[2], 1.0);
    EXPECT_NEAR(v[3], 0.0, 1.0e-12);

    EXPECT_FLOAT_EQ(v[4], 0.0);
    EXPECT_FLOAT_EQ(v[5], 1.0);
  }

  {
    mock_cell triangle;
    triangle.fill_reference_triangle(2);

    triangle_1_finite_element_0 e;
    std::vector<double> X = {0.2, 0.5};
    std::vector<double> v((X.size()/2)*3*2);

    e.evaluate_basis_derivatives_all(1, v.data(), X.data(),
                                 triangle.coordinate_dofs, 0);
    EXPECT_FLOAT_EQ(v[0], -1.0);
    EXPECT_FLOAT_EQ(v[1], -1.0);

    EXPECT_FLOAT_EQ(v[2], 1.0);
    EXPECT_NEAR(v[3], 0.0, 1.0e-12);

    EXPECT_FLOAT_EQ(v[4], 0.0);
    EXPECT_FLOAT_EQ(v[5], 1.0);
  }
}
*/

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
