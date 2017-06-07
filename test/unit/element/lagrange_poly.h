#ifndef __LAGRANGE_POLY_H
#define __LAGRANGE_POLY_H


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
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

  std::vector<std::vector<int>> poly_basis(int degree, int gdim)
  {
    std::vector<std::vector<int>> monomials;
    for (int i = 0; i <= degree ; ++i)
    {
      for (int j = 0; j <= degree ; ++j)
      {
        if ((i + j) <= degree)
        {
          monomials.push_back({i, j});
        }
      }
    }
    return monomials;
  }


  MatrixXd vandermonde(const std::vector<std::vector<int>>& monomials,
                       const boost::multi_array<double, 2>& X)
  {
    // Size matrix
    const std::size_t dim = monomials.size();
    MatrixXd A = MatrixXd::Zero(dim, dim);
    for (std::size_t i = 0; i < dim;  ++i)
    {
      // Get point
      const auto xp = X[i];

      // Fill columns for current point
      for (std::size_t j = 0; j < dim; ++j)
      {
        A(i, j) = 1.0;
        for (std::size_t d = 0; d < 2; ++d)
          A(i, j) *= std::pow(xp[d], monomials[j][d]);
      }
    }

    return A;
  }

  std::vector<double> eval(const boost::multi_array<double, 2>& points,
                           const std::vector<std::vector<int>>& monomials,
                           const std::vector<double>& coefficients)
  {
    const std::size_t gdim = 2;

    // Iterate over points
    std::vector<double> f;
    for (std::size_t i = 0; i < points.shape()[0]; ++i)
    {
      const auto x = points[i];

      // Iterate over monomial terms
      double _f = 0.0;
      for (std::size_t j = 0; j < coefficients.size(); ++j)
      {
        // Iterate over space dimensions
        double m = 1.0;
        for (std::size_t d = 0; d < gdim; ++d)
          m *= std::pow(x[d], monomials[j][d]);

        _f += coefficients[j]*m;
      }

      f.push_back(_f);
    }

    return f;
  }

  std::vector<double> evaluate_basis(const boost::multi_array<double, 2>& x,
                                     const boost::multi_array<double, 2>& vertices,
                                     int degree)
  {
    int gdim = 2;

    // Build monomials
    const std::vector<std::vector<int>> p = poly_basis(degree, gdim);
    for (auto d : p)
    {
      std::cout << "Mono: " << std::endl;
      for (auto e : d)
      {
        std::cout << "    " << e << std::endl;
      }
    }


    // Build Vandermonde matrix
    const MatrixXd A = vandermonde(p, vertices);

    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    std::cout << A.format(CleanFmt) << std::endl;

    // Factorize Vandermonde matrix
    Eigen::FullPivLU<MatrixXd> LU(A);

    std::vector<double> f;
    for (std::size_t i = 0; i < 3; ++i)
    {
      VectorXd b = VectorXd::Zero(3);
      b[i] = 1.0;

      const VectorXd c = LU.solve(b);
      for (std::size_t k = 0; k < c.rows(); ++k)
        std::cout << "C: " << c[k] << std::endl;

      std::vector<double> coeff(c.data(), c.data() + c.rows() * c.cols());
      std::vector<double> fp = eval(x, p, coeff);

      f.insert(std::end(f), std::begin(fp), std::end(fp));
    }


    return f;
  }


}

#endif
