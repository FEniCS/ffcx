#include <gtest/gtest.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <ufc.h>
#include <ufc_geometry.h>
#include <boost/multi_array.hpp>
#include "lagrange_elements.h"
#include "../../uflacs/crosslanguage/cppsupport/mock_cells.h"

#include "lagrange_poly.h"

namespace
{
  // Return list of elements
  std::vector<std::shared_ptr<ufc::finite_element>>
  get_lagrange_elements(ufc::shape shape)
  {
    std::vector<std::shared_ptr<ufc::finite_element>> elements;
    if (shape == ufc::shape::interval)
    {
      elements = {std::make_shared<interval_1_finite_element_0>(),
                  std::make_shared<interval_2_finite_element_0>(),
                  std::make_shared<interval_3_finite_element_0>(),
                  std::make_shared<interval_4_finite_element_0>()};

    }
    else if (shape == ufc::shape::triangle)
    {
      elements = {std::make_shared<triangle_1_finite_element_0>(),
                  std::make_shared<triangle_2_finite_element_0>(),
                  std::make_shared<triangle_3_finite_element_0>(),
                  std::make_shared<triangle_4_finite_element_0>()};
    }
    else if (shape == ufc::shape::tetrahedron)
    {
      elements = {std::make_shared<tetrahedron_1_finite_element_0>(),
                  std::make_shared<tetrahedron_2_finite_element_0>(),
                  std::make_shared<tetrahedron_3_finite_element_0>(),
                  std::make_shared<tetrahedron_4_finite_element_0>()};
    }
    else
      std::cerr << "Unsupported cell type" << std::endl;

    return elements;
  }


  // Tests basis evaluations on a reference cell (using
  // ufc::finite_element::evaluate_basis_reference and
  // ufc::finite_element::evaluate_basis_aa)
  void test_eval_basis_reference(const ufc::finite_element& e,
                                 const boost::multi_array<double, 2>& X)
  {
    // Number of points
    const std::size_t num_points = X.shape()[0];

    const std::size_t dim = e.space_dimension();
    const std::size_t gdim = e.geometric_dimension();
    const std::size_t degree = e.degree();
    const std::size_t ref_value_size = e.reference_value_size();

    // Get dof coordinates on reference cell
    boost::multi_array<double, 2> v(boost::extents[dim][gdim]);
    e.tabulate_reference_dof_coordinates(v.data());

    // Compute reference values
    auto f_ref = lagrange::evaluate_basis(X, v, degree, 0);

    // Compute values via FFC (reference version)
    boost::multi_array<double, 3> f(boost::extents[num_points][dim][ref_value_size]);
    std::fill_n(f.data(), f.num_elements(), 100.0);
    e.evaluate_reference_basis(f.data(), X.shape()[0], X.data());

    // Check values
    // Loop over basis functions
    for (std::size_t i = 0; i < f_ref.shape()[0]; ++i)
    {
      // Loop over points
      for (std::size_t j = 0; j < f_ref.shape()[1]; ++j)
        EXPECT_NEAR(f_ref[i][j][0], f[j][i][0], 1.0e-13);
    }

    // Test real space ufc::finite_element::evaluate_basis
    mock_cell cell;
    if (e.cell_shape() == ufc::shape::interval)
      cell.fill_reference_interval(gdim);
    else if (e.cell_shape() == ufc::shape::triangle)
      cell.fill_reference_triangle(gdim);
    else if (e.cell_shape() == ufc::shape::tetrahedron)
      cell.fill_reference_tetrahedron(gdim);
    else
      std::cerr << "Wrong cell type" << std::endl;

    // Loop over points
    std::vector<double> f_eval(dim);
    for (std::size_t j = 0; j < f_ref.shape()[1]; ++j)
    {
      std::vector<double> p(X[j].begin(), X[j].end());
      e.evaluate_basis_all(f_eval.data(), p.data(), cell.coordinate_dofs, 0);

      // Loop over basis functions
      for (std::size_t i = 0; i < f_ref.shape()[0]; ++i)
        EXPECT_NEAR(f_ref[i][j][0], f_eval[i], 1.0e-13);
    }
  }


  // Tests computes derivatives on a reference cell (using reference
  // and physical version)
  void test_reference_derivatives(const ufc::finite_element& e,
                                  const boost::multi_array<double, 2>& X)
  {
    // Number of points
    const std::size_t num_points = X.shape()[0];

    const std::size_t dim = e.space_dimension();
    const std::size_t gdim = e.geometric_dimension();
    const std::size_t degree = e.degree();
    const std::size_t ref_value_size = e.reference_value_size();

    // Get dof coordinates on reference cell
    boost::multi_array<double, 2> v(boost::extents[dim][gdim]);
    e.tabulate_reference_dof_coordinates(v.data());

    // Compute reference values
    auto f_ref = lagrange::evaluate_basis(X, v, degree, 1);

    // Compute values via FFC (reference version)
    boost::multi_array<double, 4> f(boost::extents[num_points][dim][gdim][ref_value_size]);
    e.evaluate_reference_basis_derivatives(f.data(), 1, num_points, X.data());

    // Check values
    // Loop over basis functions
    for (std::size_t i = 0; i < f_ref.shape()[0]; ++i)
    {
      // Loop over points
      for (std::size_t j = 0; j < f_ref.shape()[1]; ++j)
      {
        // Loop over derivative components
        for (std::size_t r = 0; r < f_ref.shape()[2]; ++r)
          EXPECT_NEAR(f_ref[i][j][r], f[j][i][r][0], 1.0e-13);
      }
    }

    // Test real space ufc::finite_element::evaluate_basis_derivative
    mock_cell cell;
    if (e.cell_shape() == ufc::shape::interval)
      cell.fill_reference_interval(gdim);
    else if (e.cell_shape() == ufc::shape::triangle)
      cell.fill_reference_triangle(gdim);
    else if (e.cell_shape() == ufc::shape::tetrahedron)
      cell.fill_reference_tetrahedron(gdim);
    else
      std::cerr << "Wrong cell type" << std::endl;

    // Loop over points
    std::vector<double> f_eval(dim*gdim);
    for (std::size_t p = 0; p < f_ref.shape()[1]; ++p)
    {
      std::vector<double> point(X[p].begin(), X[p].end());
      e.evaluate_basis_derivatives_all(1, f_eval.data(), point.data(),
                                       cell.coordinate_dofs, 0);

      // Loop over basis functions
      for (std::size_t d = 0; d < f_ref.shape()[0]; ++d)
      {
        // Loop over derivative components
        for (std::size_t r = 0; r < f_ref.shape()[2]; ++r)
          EXPECT_NEAR(f_ref[d][p][r], f_eval[gdim*d + r], 1.0e-12);
      }
    }
  }
}

/*
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
*/



TEST(FiniteElementScalarLagrange, eval_basis)
{
  // Interval elements
  {
    // Points at which to evaluate basis
    boost::multi_array<double, 2> X(boost::extents[4][1]);
    X[0][0] = 0.5;  X[1][0] = 1.0; X[2][0] = 0.5; X[3][0] = 0.5;

    // Lists of elements and test
    auto elements =  get_lagrange_elements(ufc::shape::interval);
    for (auto e : elements)
    {
      assert(e);
      test_eval_basis_reference(*e, X);
    }
  }

  // Triangles
  {
    // Points at which to evaluate basis
    boost::multi_array<double, 2> X(boost::extents[4][2]);
    X[0][0] = 0.5;  X[0][1] = 0.0;
    X[1][0] = 1.0;  X[1][1] = 0.0;
    X[2][0] = 0.5;  X[2][1] = 0.0;
    X[3][0] = 0.5;  X[3][1] = 0.5;

    // Get lists of elements and test
    auto elements =  get_lagrange_elements(ufc::shape::triangle);
    for (auto e : elements)
    {
      assert(e);
      test_eval_basis_reference(*e, X);
    }
  }

  // Tetrahedra
  {
    // Points at which to evaluate basis
    boost::multi_array<double, 2> X(boost::extents[4][3]);
    X[0][0] = 0.5;  X[0][1] = 0.0; X[0][2] = 0.0;
    X[1][0] = 1.0;  X[1][1] = 0.0; X[0][2] = 1.0;
    X[2][0] = 0.5;  X[2][1] = 0.0; X[0][2] = 0.0;
    X[3][0] = 0.5;  X[3][1] = 0.5; X[0][2] = 0.2;

    // Get lists of elements and test
    auto elements =  get_lagrange_elements(ufc::shape::tetrahedron);
    for (auto e : elements)
    {
      assert(e);
      test_eval_basis_reference(*e, X);
    }
  }
}


TEST(FiniteElementScalarLagrange, eval_basis_d)
{
  // Interval elements
  {
    // Points at which to evaluate basis
    boost::multi_array<double, 2> X(boost::extents[4][1]);
    X[0][0] = 0.5;  X[1][0] = 1.0; X[2][0] = 0.5; X[3][0] = 0.5;

    // Get lists of elements and test
    auto elements =  get_lagrange_elements(ufc::shape::interval);
    for (auto e : elements)
    {
      assert(e);
      test_reference_derivatives(*e, X);
    }
  }

  // Triangles
  {
    // Points at which to evaluate basis
    boost::multi_array<double, 2> X(boost::extents[4][2]);
    X[0][0] = 0.5;  X[0][1] = 0.0;
    X[1][0] = 1.0;  X[1][1] = 0.0;
    X[2][0] = 0.5;  X[2][1] = 0.0;
    X[3][0] = 0.5;  X[3][1] = 0.5;

    // Get lists of elements and test
    auto elements =  get_lagrange_elements(ufc::shape::triangle);
    for (auto e : elements)
    {
      assert(e);
      test_reference_derivatives(*e, X);
    }
  }

  // Tetrahedra
  {
    // Points at which to evaluate basis
    boost::multi_array<double, 2> X(boost::extents[4][3]);
    X[0][0] = 0.5;  X[0][1] = 0.0; X[0][2] = 0.0;
    X[1][0] = 1.0;  X[1][1] = 0.0; X[0][2] = 1.0;
    X[2][0] = 0.5;  X[2][1] = 0.0; X[0][2] = 0.0;
    X[3][0] = 0.5;  X[3][1] = 0.5; X[0][2] = 0.2;

    // Get lists of elements and test
    auto elements =  get_lagrange_elements(ufc::shape::tetrahedron);
    for (auto e : elements)
    {
      assert(e);
      test_reference_derivatives(*e, X);
    }
  }
}

TEST(FiniteElementScalarLagrange, cell_shape)
{
  // Shapes
  std::set<ufc::shape> shapes = {ufc::shape::interval,
                                 ufc::shape::triangle,
                                 ufc::shape::tetrahedron};

  for (auto shape : shapes)
  {
    auto elements = get_lagrange_elements(shape);
    for (auto e : elements)
    {
      assert(e);
      EXPECT_EQ(e->cell_shape(), shape);
    }
  }
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
