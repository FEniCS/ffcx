// This is utility code for UFC (Unified Form-assembly Code) v 2.0.5.
// This code is released into the public domain.
// Modified for UFLACS test suite 2012.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2006-2012.

#ifndef __UFC_CELLS_H
#define __UFC_CELLS_H

#include "ufc.h"
#include <stdexcept>

namespace uflacs
{
    unsigned shape_to_dimension(ufc::shape s)
    {
        switch (s)
        {
        case ufc::interval:
            return 1;
        case ufc::triangle:
            return 2;
        case ufc::tetrahedron:
            return 3;
        case ufc::quadrilateral:
            return 2;
        case ufc::hexahedron:
            return 3;
        }
        return 0;
    }

  class mock_cell
  {
  public:

    /// Constructor
    cell(ufc::shape s):
        cell_shape(s),
        topological_dimension(shape_to_dimension(s)),
        geometric_dimension(topological_dimension),
        local_facet(-1)
    {
    }

    /// Shape of the cell
    shape cell_shape;

    /// Topological dimension of the mesh
    unsigned int topological_dimension;

    /// Geometric dimension of the mesh
    unsigned int geometric_dimension;

    /// Array of coordinates for the vertices of the cell, large enough for all valid cells
    double coordinates[8][3];

    /// Local facet index
    int local_facet;
  };


    void allocate_cell_arrays(ufc::cell & c, int *num_entities)
    {
        // Fill global indices like we had a single-cell mesh.
        c.entity_indices = new unsigned int*[c.topological_dimension+1];
        for(unsigned int i = 0; i <= c.topological_dimension; ++i)
        {
            c.entity_indices[i] = new unsigned int[num_entities[i]];
            for(unsigned int j = 0; j < num_entities[i]; ++j)
            {
                c.entity_indices[i][j] = j;
            }
        }

        // Allocate an empty array of vertex coordinates.
        c.coordinates = new double*[num_entities[0]];
        for(unsigned int i = 0; i < num_entities[0]; ++i)
        {
            c.coordinates[i] = new double[c.geometric_dimension];
            for(unsigned int j = 0; j < c.geometric_dimension; ++j)
            {
                c.coordinates[i][j] = 0.0;
            }
        }
    }

    void make_reference_cell(ufc::shape s, ufc::cell & c)
    {
        switch(s)
        {
        case ufc::interval:
            make_reference_interval(c);
            break;
        case ufc::triangle:
            make_reference_triangle(c);
            break;
        case ufc::quadrilateral:
            make_reference_quadrilateral(c);
            break;
        case ufc::tetrahedron:
            make_reference_tetrahedron(c);
            break;
        case ufc::hexahedron:
            make_reference_hexahedron(c);
            break;
        default:
            throw std::runtime_error("Invalid shape.");
        }
    }

    void make_reference_interval(ufc::cell & c)
    {
        // Configure dimensions
        c.cell_shape = ufc::interval;
        c.topological_dimension = 1;
        c.geometric_dimension = c.topological_dimension;
        int num_entities[2] = { 2, 1 };

        // Allocate arrays consistent with dimensions
        allocate_cell_arrays(c, num_entities);

        // Fill in coordinates for reference cell
        coordinates[0][0] = 0.0;
        coordinates[1][0] = 1.0;
    }

    void make_reference_triangle(ufc::cell & c)
    {
        // Configure dimensions
        c.cell_shape = ufc::triangle;
        c.topological_dimension = 2;
        c.geometric_dimension = c.topological_dimension;
        int num_entities[3] = { 3, 3, 1 };

        // Allocate arrays consistent with dimensions
        allocate_cell_arrays(c, num_entities);

        // Fill in coordinates for reference cell
        coordinates[0][0] = 0.0;
        coordinates[0][1] = 0.0;

        coordinates[1][0] = 1.0;
        coordinates[1][1] = 0.0;

        coordinates[2][0] = 0.0;
        coordinates[2][1] = 1.0;
    }

    void make_reference_tetrahedron(ufc::cell & c)
    {
        // Configure dimensions
        c.cell_shape = ufc::tetrahedron;
        c.topological_dimension = 3;
        c.geometric_dimension = c.topological_dimension;
        int num_entities[4] = { 4, 6, 4, 1 };

        // Allocate arrays consistent with dimensions
        allocate_cell_arrays(c, num_entities);

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

    void make_reference_quadrilateral(ufc::cell & c)
    {
        // Configure dimensions
        c.cell_shape = ufc::quadrilateral;
        c.topological_dimension = 2;
        c.geometric_dimension = c.topological_dimension;
        int num_entities[3] = { 4, 4, 1 };

        // Allocate arrays consistent with dimensions
        allocate_cell_arrays(c, num_entities);

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

    void make_reference_hexahedron(ufc::cell & c)
    {
        // Configure dimensions
        c.cell_shape = ufc::hexahedron;
        c.topological_dimension = 3;
        c.geometric_dimension = c.topological_dimension;
        int num_entities[4] = { 8, 12, 6, 1 };

        // Allocate arrays consistent with dimensions
        allocate_cell_arrays(c, num_entities);

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

    void scale_cell_coordinates(Cell & c, double *factor)
    {
        for (int i = 0; i < c.cell_)
    }

    void scale_cell_coordinates(cell & c, double factor)
    {
        double vfactor[3] = { factor, factor, factor };
        scale_cell_coordinates(c, vfactor);
    }

    /// Convenience subclass of ufc::cell for constructing a cell from data
    class Cell: public ufc::cell
    {
    public:

        Cell(ufc::shape)
        {
            make_reference_cell(self);
        }

        Cell(unsigned int top, unsigned int geo,
             std::vector< std::vector<double> > coords,
             std::vector< unsigned int> num_ents):
            ufc::cell(), num_entities(num_ents)
        {
            topological_dimension = top;
            geometric_dimension   = geo;
            num_entities[0] = coords.size();

            // Fill global indices
            entity_indices = new unsigned int*[topological_dimension+1];
            for(unsigned int i = 0; i <= topological_dimension; ++i)
            {
                entity_indices[i] = new unsigned int[num_entities[i]];
                for(unsigned int j = 0; j < num_entities[i]; ++j)
                {
                    entity_indices[i][j] = j;
                }
            }

            for(unsigned int i = 0; i < num_ents.size(); ++i)
            {
                num_entities[i] = num_ents[i];
            }

            // Allocate an empty array of vertex coordinates.
            coordinates = new double*[coords.size()];
            for(unsigned int i = 0; i < coords.size(); ++i)
            {
                coordinates[i] = new double[geometric_dimension];
                for(unsigned int j = 0; j < geometric_dimension; ++j)
                {
                    coordinates[i][j] = coords[i][j];
                }
            }

        }

        /// Destructor
        virtual ~Cell()
        {
            for(unsigned int i = 0; i <= topological_dimension; ++i)
            {
                delete [] entity_indices[i];
            }
            delete [] entity_indices;

            for(unsigned int i = 0; i < num_entities[0]; ++i)
            {
                delete [] coordinates[i];
            }
            delete [] coordinates;
        }

        /// The number of entities of a particular dimension
        std::vector<unsigned int> num_entities;
    };

}

#endif
