#ifndef __LAGRANGE_POLY_H
#define __LAGRANGE_POLY_H


#include <cassert>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <ufc.h>
#include <ufc_geometry.h>
#include <boost/multi_array.hpp>
#include "lagrange_elements.h"
#include "../../uflacs/crosslanguage/cppsupport/mock_cells.h"



namespace lagrange
{
  // This collection of function computes Lagrange polynomial on
  // simplices via the Vandermonde matrix. It is used to test UFC
  // code, and the approach is not advocated for production runs. It
  // is used to test the FFC-genereted basis evaluation code.

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;


  // Compute monomial exponents for complete polynomial
  boost::multi_array<unsigned int, 2> poly_basis(int degree, int gdim)
  {
    std::vector<std::vector<int>> m;
    switch (gdim)
    {
    case 1:
      for (int i = 0; i <= degree ; ++i)
        m.push_back({i});
      break;
    case 2:
      for (int i = 0; i <= degree ; ++i)
      {
        for (int j = 0; j <= degree ; ++j)
        {
          if ((i + j) <= degree)
            m.push_back({i, j});
        }
      }
      break;
    case 3:
      for (int i = 0; i <= degree ; ++i)
      {
        for (int j = 0; j <= degree ; ++j)
        {
          for (int k = 0; k <= degree ; ++k)
          {
            if ((i + j + k) <= degree)
              m.push_back({i, j, k});
          }
        }
      }
      break;
    default:
      std::cerr << "poly_basis support 1D, 2D and 3D only." << std::endl;
    }

    boost::multi_array<unsigned int, 2> monomials(boost::extents[m.size()][gdim]);
    for (std::size_t i = 0; i < m.size(); ++i)
      for (int j = 0; j < gdim; ++j)
        monomials[i][j] = m[i][j];

    return monomials;
  }


  // Build Vandermonde matrix
  MatrixXd vandermonde(const boost::multi_array<unsigned int, 2>& monomials,
                       const boost::multi_array<double, 2>& X)
  {
    const std::size_t dim = monomials.shape()[0];
    const std::size_t gdim = X.shape()[1];
    assert(X.shape()[0] == dim);
    assert(monomials.shape()[1] == gdim);

    // Create  matrix and build
    MatrixXd A = MatrixXd::Zero(dim, dim);
    for (std::size_t i = 0; i < dim;  ++i)
    {
      // Get point for this row
      const auto xp = X[i];
      //std::cout << "Point (a) : " << xp[0] << ", " << xp[1] << std::endl;
      //std::cout << "Point (b) : " << X[i][0] << ", " << X[i][1] << std::endl;

      // Fill columns for current point
      for (std::size_t j = 0; j < dim; ++j)
      {
        A(i, j) = 1.0;
        const auto& m = monomials[j];
        for (std::size_t d = 0; d < gdim; ++d)
          A(i, j) *= std::pow(xp[d], m[d]);
      }
    }

    return A;
  }

  // Evaluate a polynomial at collection of points
  std::vector<double> eval(const boost::multi_array<double, 2>& points,
                           const boost::multi_array<unsigned int, 2>& monomials,
                           const std::vector<double>& coefficients)
  {
    // Geometric dimension
    const std::size_t gdim = points.shape()[1];
    assert(!monomials.empty());
    assert(monomials[0].size() == gdim);

    // Vector to hold values at each point
    std::vector<double> f;

    // Iterate over points and evaluate polynomial at each point
    for (std::size_t i = 0; i < points.shape()[0]; ++i)
    {
      const auto x = points[i];

      // Iterate over monomial terms
      double _f = 0.0;
      for (std::size_t j = 0; j < coefficients.size(); ++j)
      {
        const auto& exponents = monomials[j];

        // Iterate over space dimensions
        double m = 1.0;
        for (std::size_t d = 0; d < gdim; ++d)
          m *= std::pow(x[d], exponents[d]);

        _f += coefficients[j]*m;
      }

      f.push_back(_f);
    }

    return f;
  }

  // Evaluates Lagrange (scalar) basis functions at points on affine
  // cells, returning f[basis_index][point]
  boost::multi_array<double, 2> evaluate_basis(const boost::multi_array<double, 2>& x,
                                               const boost::multi_array<double, 2>& vertices,
                                               int degree)
  {
    // Geometric dim
    const std::size_t gdim = x.shape()[1];
    assert(gdim == vertices.shape()[1]);

    // Number of points
    const std::size_t num_points = x.shape()[0];

    // Build monomials
    const boost::multi_array<unsigned int, 2> p = poly_basis(degree, gdim);
    /*
    for (auto d : p)
    {
      std::cout << "Mono: " << std::endl;
      for (auto e : d)
      {
        std::cout << "    " << e << std::endl;
      }
    }
    */

    // Dimension of polynomial space
    const std::size_t dim = p.shape()[0];

    // Build Vandermonde matrix
    const MatrixXd A = vandermonde(p, vertices);
    //Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    //std::cout << A.format(CleanFmt) << std::endl;

    // Factorize Vandermonde matrix
    Eigen::FullPivLU<MatrixXd> LU(A);

    // Iterate over basis functions and evaluate at points
    boost::multi_array<double, 2> f(boost::extents[dim][num_points]);
    for (std::size_t i = 0; i < dim; ++i)
    {
      VectorXd b = VectorXd::Zero(dim);
      b[i] = 1.0;

      // Solve and copy coefficients into a vector
      const VectorXd c = LU.solve(b);
      std::vector<double> coeff(c.data(), c.data() + c.rows());
      //for (std::size_t k = 0; k < c.rows(); ++k)
      //  std::cout << "C: " << c[k] << std::endl;

      // Evaluate basis function i at each point
      std::vector<double> fp = eval(x, p, coeff);

      // Copy result in top return array
      std::copy(std::begin(fp), std::end(fp), f[i].begin());

      //for (auto _f : fp)
      //  std::cout << "   Basis function eval: " << _f << std::endl;

      //f.insert(std::end(f), std::begin(fp), std::end(fp));
    }

    return f;
  }


}

#endif
