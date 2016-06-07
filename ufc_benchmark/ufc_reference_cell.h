// This is utility code for UFC (Unified Form-assembly Code).
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2006-2015.

#ifndef __UFC_REFERENCE_CELL_H
#define __UFC_REFERENCE_CELL_H

#include "ufc.h"
#include <cstddef>
#include <stdexcept>

namespace ufc
{

    /// Description of a reference cell, for debugging and testing UFC code.
    class reference_cell: public ufc::cell
    {
    public:

        /// Constructor
        reference_cell(ufc::shape s)
        {
            cell_shape = s;

            num_entities[0] = 0;
            num_entities[1] = 0;
            num_entities[2] = 0;
            num_entities[3] = 0;

            // Get topological dimension and number of entities in a cell of this type.
            switch(s)
            {
            case interval:      topological_dimension = 1; num_entities[0] = 2; num_entities[1] = 1;  break;
            case triangle:      topological_dimension = 2; num_entities[0] = 3; num_entities[1] = 3;  num_entities[2] = 1; break;
            case quadrilateral: topological_dimension = 2; num_entities[0] = 4; num_entities[1] = 4;  num_entities[2] = 1; break;
            case tetrahedron:   topological_dimension = 3; num_entities[0] = 4; num_entities[1] = 6;  num_entities[2] = 4; num_entities[3] = 1; break;
            case hexahedron:    topological_dimension = 3; num_entities[0] = 8; num_entities[1] = 12; num_entities[2] = 6; num_entities[3] = 1; break;
            default: throw std::runtime_error("Invalid shape.");
            }

            // Assume same geometric dimension.
            geometric_dimension = topological_dimension;

            // Fill global indices like we had a single-cell mesh.
            entity_indices = new std::size_t*[topological_dimension+1];
            for(std::size_t i = 0; i <= topological_dimension; i++)
            {
                entity_indices[i] = new std::size_t[num_entities[i]];
                for(std::size_t j = 0; j < num_entities[i]; j++)
                {
                    entity_indices[i][j] = j;
                }
            }

            // Allocate an empty array of vertex coordinates.
            coordinates = new double*[num_entities[0]];
            for(std::size_t i = 0; i < num_entities[0]; i++)
            {
                coordinates[i] = new double[geometric_dimension];
                for(std::size_t j = 0; j < geometric_dimension; j++)
                {
                    coordinates[i][j] = 0.0;
                }
            }

            // Fill coordinates with reference cell definition.
            switch(s)
            {
            case interval:
                coordinates[0][0] = 0.0;
                coordinates[1][0] = 1.0;
                break;

            case triangle:
                coordinates[0][0] = 0.0;
                coordinates[0][1] = 0.0;

                coordinates[1][0] = 1.0;
                coordinates[1][1] = 0.0;

                coordinates[2][0] = 0.0;
                coordinates[2][1] = 1.0;
                break;

            case quadrilateral:
                coordinates[0][0] = 0.0;
                coordinates[0][1] = 0.0;

                coordinates[1][0] = 1.0;
                coordinates[1][1] = 0.0;

                coordinates[2][0] = 1.0;
                coordinates[2][1] = 1.0;

                coordinates[3][0] = 0.0;
                coordinates[3][1] = 1.0;
                break;

            case tetrahedron:
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
                break;

            case hexahedron:
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
                break;
            }
        }

        /// Destructor
        virtual ~reference_cell()
        {
            for(std::size_t i = 0; i <= topological_dimension; i++)
            {
                delete [] entity_indices[i];
            }
            delete [] entity_indices;

            for(std::size_t i = 0; i < num_entities[0]; i++)
            {
                delete [] coordinates[i];
            }
            delete [] coordinates;
        }

        /// The number of entities of a particular dimension
        std::size_t num_entities[4];

    };

    /// Description of a reference cell, for debugging and testing UFC code.
    class Cell: public ufc::cell
    {
    public:

        /// Constructor
        Cell(std::size_t top, std::size_t geo, std::vector< std::vector<double> > coords, std::vector< std::size_t> num_ents): ufc::cell(), num_entities(num_ents)
        {
            topological_dimension = top;
            geometric_dimension   = geo;
            num_entities[0] = coords.size();

            // Fill global indices
//            entity_indices = new std::size_t*[topological_dimension+1];
//            for(std::size_t i = 0; i <= topological_dimension; i++)
//            {
//                entity_indices[i] = new std::size_t[num_entities[i]];
//                for(std::size_t j = 0; j < num_entities[i]; j++)
//                {
//                    entity_indices[i][j] = j;
//                }
//            }

            for(std::size_t i = 0; i < num_ents.size(); i++)
              num_entities[i] = num_ents[i];

            // Allocate an empty array of vertex coordinates.
            coordinates = new double*[coords.size()];
            for(std::size_t i = 0; i < coords.size(); i++)
            {
                coordinates[i] = new double[geometric_dimension];
                for(std::size_t j = 0; j < geometric_dimension; j++)
                {
                    coordinates[i][j] = coords[i][j];
                }
            }

        }

        /// Destructor
        virtual ~Cell()
        {
//            for(std::size_t i = 0; i <= topological_dimension; i++)
//            {
//                delete [] entity_indices[i];
//            }
//            delete [] entity_indices;

            for(std::size_t i = 0; i < num_entities[0]; i++)
            {
                delete [] coordinates[i];
            }
            delete [] coordinates;
        }

        /// The number of entities of a particular dimension
        std::vector<std::size_t> num_entities;
    };


    /// Consistent data for a mesh consisting of a single reference cell, for debugging and testing UFC code.
    class reference_mesh: public ufc::mesh
    {
    public:

        /// Constructor
        reference_mesh(ufc::shape s):
            c(s)
        {
            topological_dimension = c.topological_dimension;
            geometric_dimension   = c.geometric_dimension;

            // Set global number of entities of each topological dimension to that of a single cell.
            num_entities = new std::size_t[topological_dimension+1];
            for(std::size_t i = 0; i <= topological_dimension; i++)
            {
                num_entities[i] = c.num_entities[i];
            }
        }

        /// Destructor
        virtual ~reference_mesh()
        {
            delete [] num_entities;
        }

        /// A reference cell, the only cell in this mesh.
        reference_cell c;

    };

    /// Consistent data for a mesh consisting of a single reference cell, for debugging and testing UFC code.
    class Mesh: public ufc::mesh
    {
    public:

        /// Constructor
        Mesh(std::size_t top, std::size_t geo, std::vector<std::size_t> ents)//: ufc::mesh()
        {
            topological_dimension = top;
            geometric_dimension   = geo;

            // Set global number of entities of each topological dimension to that of a single cell.
            num_entities = new std::size_t[topological_dimension+1];
            for(std::size_t i = 0; i <= topological_dimension; i++)
            {
                num_entities[i] = ents[i];
            }
        }

        /// Destructor
        virtual ~Mesh()
        {
            delete [] num_entities;
        }
    };

}

#endif
