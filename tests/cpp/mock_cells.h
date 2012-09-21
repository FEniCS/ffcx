// This is utility code for UFLACS unit tests.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2006-2012.

#ifndef __MOCK_CELLS_H
#define __MOCK_CELLS_H

#include "ufc.h"

namespace uflacs
{
    struct mock_cell
    {
        mock_cell(ufc::shape cell_shape, size_t geometric_dimension,
                  size_t topological_dimension, size_t num_vertices):
          cell_shape(cell_shape),
          geometric_dimension(geometric_dimension),
          topological_dimension(topological_dimension),
          num_vertices(num_vertices),
          local_facet(0)
        {
          for (int i=0; i<12; ++i)
            for (int j=0; j<3; ++j)
              coordinates[i][j] = 0.0;
        }

        ufc::shape cell_shape;
        size_t geometric_dimension;
        size_t topological_dimension;
        size_t num_vertices; // not part of ufc::cell but needed for cell mapping code
        int local_facet;
        double coordinates[12][3];
    };

    struct mock_interval: public mock_cell
    {
        mock_interval():
            mock_cell(ufc::interval, 1, 1, 2)
        {
            // Fill in coordinates for reference cell
            coordinates[0][0] = 0.0;
            coordinates[1][0] = 1.0;
        }
    };

    struct mock_triangle: public mock_cell
    {
        mock_triangle():
            mock_cell(ufc::triangle, 2, 2, 3)
        {
            // Fill in coordinates for reference cell
            coordinates[0][0] = 0.0;
            coordinates[0][1] = 0.0;

            coordinates[1][0] = 1.0;
            coordinates[1][1] = 0.0;

            coordinates[2][0] = 0.0;
            coordinates[2][1] = 1.0;
        }
    };

    struct mock_tetrahedron: public mock_cell
    {
        mock_tetrahedron():
            mock_cell(ufc::tetrahedron, 3, 3, 4)
        {
            // Fill in coordinates for reference cell
            coordinates[0][0] = 0.0;
            coordinates[0][1] = 0.0;
            coordinates[0][2] = 0.0;

            coordinates[1][0] = 1.0;
            coordinates[1][1] = 0.0;
            coordinates[1][2] = 0.0;

            coordinates[2][0] = 0.0;
            coordinates[2][1] = 1.0;
            coordinates[2][2] = 0.0;

            coordinates[3][0] = 0.0;
            coordinates[3][1] = 0.0;
            coordinates[3][2] = 1.0;
        }
    };

    struct mock_quadrilateral: public mock_cell
    {
        mock_quadrilateral():
            mock_cell(ufc::quadrilateral, 2, 2, 4)
        {
            // Fill in coordinates for reference cell
            coordinates[0][0] = 0.0;
            coordinates[0][1] = 0.0;

            coordinates[1][0] = 1.0;
            coordinates[1][1] = 0.0;

            coordinates[2][0] = 1.0;
            coordinates[2][1] = 1.0;

            coordinates[3][0] = 0.0;
            coordinates[3][1] = 1.0;
        }
    };

    struct mock_hexahedron: public mock_cell
    {
        mock_hexahedron():
            mock_cell(ufc::hexahedron, 3, 3, 8)
        {
            // Fill in coordinates for reference cell
            coordinates[0][0] = 0.0;
            coordinates[0][1] = 0.0;
            coordinates[0][2] = 0.0;

            coordinates[1][0] = 1.0;
            coordinates[1][1] = 0.0;
            coordinates[1][2] = 0.0;

            coordinates[2][0] = 1.0;
            coordinates[2][1] = 1.0;
            coordinates[2][2] = 0.0;

            coordinates[3][0] = 0.0;
            coordinates[3][1] = 1.0;
            coordinates[3][2] = 0.0;

            coordinates[4][0] = 0.0;
            coordinates[4][1] = 0.0;
            coordinates[4][2] = 1.0;

            coordinates[5][0] = 1.0;
            coordinates[5][1] = 0.0;
            coordinates[5][2] = 1.0;

            coordinates[6][0] = 1.0;
            coordinates[6][1] = 1.0;
            coordinates[6][2] = 1.0;

            coordinates[7][0] = 0.0;
            coordinates[7][1] = 1.0;
            coordinates[7][2] = 1.0;
        }
    };

    void scale_cell(mock_cell & c, double factor)
    {
        for (size_t i=0; i<c.num_vertices; ++i)
        {
            for (size_t j=0; j<c.geometric_dimension; ++j)
            {
                c.coordinates[i][j] *= factor;
            }
        }
    }

    void scale_cell(mock_cell & c, const double * factors)
    {
        for (size_t i=0; i<c.num_vertices; ++i)
        {
            for (size_t j=0; j<c.geometric_dimension; ++j)
            {
                c.coordinates[i][j] *= factors[j];
            }
        }
    }

    void translate_cell(mock_cell & c, const double * x)
    {
        for (size_t i=0; i<c.num_vertices; ++i)
        {
            for (size_t j=0; j<c.geometric_dimension; ++j)
            {
                c.coordinates[i][j] += x[j];
            }
        }
    }

    void map_cell(mock_cell & c, const double * G, const double * x0)
    {
        for (size_t i=0; i<c.num_vertices; ++i)
        {
            double result[3];
            for (size_t j=0; j<c.geometric_dimension; ++j)
            {
                result[j] = x0[j];
                for (size_t k=0; k<c.geometric_dimension; ++k)
                {
                     result[j] += G[c.geometric_dimension*j + k] * c.coordinates[i][k];
                }
            }
            for (size_t j=0; j<c.geometric_dimension; ++j)
            {
                c.coordinates[i][j] = result[j];
            }
        }
    }
}

#endif
