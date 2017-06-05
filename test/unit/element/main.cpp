
#include <gtest/gtest.h>


#include <ufc.h>
#include <ufc_geometry.h>
#include "lagrange_elements.h"
//#include "mock_cells.h"
//#include "debugging.h"


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

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
