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
        // These coordinates are all that generated code should care about, the rest is to support the tests
        double vertex_coordinates[8*3]; // Worst case hexahedron: 8 vertices in 3d

        // Dimensions needed for generic cell transformations
        size_t geometric_dimension;
        size_t topological_dimension;
        size_t num_vertices;

        // Constructor just clears everything to zero
        mock_cell()
        {
            init_dimensions(0,0,0);
        }

        // Utility initialization function
        void init_dimensions(size_t geometric_dimension, size_t topological_dimension, size_t num_vertices)
        {
            this->geometric_dimension = geometric_dimension;
            this->topological_dimension = topological_dimension;
            this->num_vertices = num_vertices;
            for (int i=0; i<sizeof(vertex_coordinates)/sizeof(vertex_coordinates[0]); ++i)
                vertex_coordinates[i] = 0.0;
        }

        void fill_reference_interval(size_t geometric_dimension)
        {
            init_dimensions(geometric_dimension, 1, 2);

            // Fill in vertex_coordinates for reference cell
            vertex_coordinates[0*geometric_dimension + 0] = 0.0;
            vertex_coordinates[1*geometric_dimension + 0] = 1.0;
        }

        void fill_reference_triangle(size_t geometric_dimension)
        {
            init_dimensions(geometric_dimension, 2, 3);

            // Fill in vertex_coordinates for reference cell
            vertex_coordinates[0*geometric_dimension + 0] = 0.0;
            vertex_coordinates[0*geometric_dimension + 1] = 0.0;

            vertex_coordinates[1*geometric_dimension + 0] = 1.0;
            vertex_coordinates[1*geometric_dimension + 1] = 0.0;

            vertex_coordinates[2*geometric_dimension + 0] = 0.0;
            vertex_coordinates[2*geometric_dimension + 1] = 1.0;
        }

        void fill_reference_tetrahedron(size_t geometric_dimension)
        {
            init_dimensions(geometric_dimension, 3, 4);

            // Fill in vertex_coordinates for reference cell
            vertex_coordinates[0*geometric_dimension + 0] = 0.0;
            vertex_coordinates[0*geometric_dimension + 1] = 0.0;
            vertex_coordinates[0*geometric_dimension + 2] = 0.0;

            vertex_coordinates[1*geometric_dimension + 0] = 1.0;
            vertex_coordinates[1*geometric_dimension + 1] = 0.0;
            vertex_coordinates[1*geometric_dimension + 2] = 0.0;

            vertex_coordinates[2*geometric_dimension + 0] = 0.0;
            vertex_coordinates[2*geometric_dimension + 1] = 1.0;
            vertex_coordinates[2*geometric_dimension + 2] = 0.0;

            vertex_coordinates[3*geometric_dimension + 0] = 0.0;
            vertex_coordinates[3*geometric_dimension + 1] = 0.0;
            vertex_coordinates[3*geometric_dimension + 2] = 1.0;
        }

        void fill_reference_quadrilateral(size_t geometric_dimension)
        {
            init_dimensions(geometric_dimension, 2, 4);

            // Fill in vertex_coordinates for reference cell
            vertex_coordinates[0*geometric_dimension + 0] = 0.0;
            vertex_coordinates[0*geometric_dimension + 1] = 0.0;

            vertex_coordinates[1*geometric_dimension + 0] = 1.0;
            vertex_coordinates[1*geometric_dimension + 1] = 0.0;

            vertex_coordinates[2*geometric_dimension + 0] = 1.0;
            vertex_coordinates[2*geometric_dimension + 1] = 1.0;

            vertex_coordinates[3*geometric_dimension + 0] = 0.0;
            vertex_coordinates[3*geometric_dimension + 1] = 1.0;
        }

        void fill_reference_hexahedron(size_t geometric_dimension)
        {
            init_dimensions(geometric_dimension, 3, 8);

            // Fill in vertex_coordinates for reference cell
            vertex_coordinates[0*geometric_dimension + 0] = 0.0;
            vertex_coordinates[0*geometric_dimension + 1] = 0.0;
            vertex_coordinates[0*geometric_dimension + 2] = 0.0;

            vertex_coordinates[1*geometric_dimension + 0] = 1.0;
            vertex_coordinates[1*geometric_dimension + 1] = 0.0;
            vertex_coordinates[1*geometric_dimension + 2] = 0.0;

            vertex_coordinates[2*geometric_dimension + 0] = 1.0;
            vertex_coordinates[2*geometric_dimension + 1] = 1.0;
            vertex_coordinates[2*geometric_dimension + 2] = 0.0;

            vertex_coordinates[3*geometric_dimension + 0] = 0.0;
            vertex_coordinates[3*geometric_dimension + 1] = 1.0;
            vertex_coordinates[3*geometric_dimension + 2] = 0.0;

            vertex_coordinates[4*geometric_dimension + 0] = 0.0;
            vertex_coordinates[4*geometric_dimension + 1] = 0.0;
            vertex_coordinates[4*geometric_dimension + 2] = 1.0;

            vertex_coordinates[5*geometric_dimension + 0] = 1.0;
            vertex_coordinates[5*geometric_dimension + 1] = 0.0;
            vertex_coordinates[5*geometric_dimension + 2] = 1.0;

            vertex_coordinates[6*geometric_dimension + 0] = 1.0;
            vertex_coordinates[6*geometric_dimension + 1] = 1.0;
            vertex_coordinates[6*geometric_dimension + 2] = 1.0;

            vertex_coordinates[7*geometric_dimension + 0] = 0.0;
            vertex_coordinates[7*geometric_dimension + 1] = 1.0;
            vertex_coordinates[7*geometric_dimension + 2] = 1.0;
        }

        // Scale cell coordinates
        void scale(double factor)
        {
            for (size_t i=0; i<num_vertices; ++i)
            {
                for (size_t j=0; j<geometric_dimension; ++j)
                {
                    vertex_coordinates[i*geometric_dimension + j] *= factor;
                }
            }
        }

        // Scale cell coordinates differently in each geometric dimension
        void scale(const double * factors)
        {
            for (size_t i=0; i<num_vertices; ++i)
            {
                for (size_t j=0; j<geometric_dimension; ++j)
                {
                    vertex_coordinates[i*geometric_dimension + j] *= factors[j];
                }
            }
        }

        // Translate cell coordinates
        void translate(const double * x)
        {
            for (size_t i=0; i<num_vertices; ++i)
            {
                for (size_t j=0; j<geometric_dimension; ++j)
                {
                    vertex_coordinates[i*geometric_dimension + j] += x[j];
                }
            }
        }

        // Apply affine mapping to cell coordinates
        void affine_map(const double * A, const double * x0)
        {
            for (size_t i=0; i<num_vertices; ++i)
            {
                double result[3];
                for (size_t j=0; j<geometric_dimension; ++j)
                {
                    result[j] = x0[j];
                    for (size_t k=0; k<geometric_dimension; ++k)
                    {
                         result[j] += A[geometric_dimension*j + k] * vertex_coordinates[i*geometric_dimension + k];
                    }
                }
                for (size_t j=0; j<geometric_dimension; ++j)
                {
                    vertex_coordinates[i*geometric_dimension + j] = result[j];
                }
            }
        }
    };

}

#endif
