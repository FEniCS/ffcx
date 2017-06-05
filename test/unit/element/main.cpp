#include <gtest/gtest.h>

#include <ufc.h>
#include <ufc_geometry.h>
#include "lagrange_elements.h"
//#include "mock_cells.h"
//#include "debugging.h"

#include "../../uflacs/crosslanguage/cppsupport/mock_cells.h"

//#include "test_ufc_integral_types.h"


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


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
