//
// This code complies with UFC version 1.0, and is generated with SyFi version 0.4.0.
//
// http://www.fenics.org/syfi/
// http://www.fenics.org/ufc/
//


#ifndef __dof_map_Lagrange_1_2D_H
#define __dof_map_Lagrange_1_2D_H

#include <stdexcept>
#include <math.h>
#include <ufc.h>
#include <pycc/Functions/Ptv.h>
#include <pycc/Functions/Ptv_tools.h>
#include <pycc/Functions/Dof_Ptv.h>
#include <pycc/Functions/OrderedPtvSet.h>
#include <pycc/Functions/Dof_OrderedPtvSet.h>



namespace pycc
{

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class dof_map_Lagrange_1_2D: public ufc::dof_map
{  public:
    pycc::Dof_Ptv dof;
    unsigned int num_elements;
    unsigned int * loc2glob;

public:

  /// Constructor
  dof_map_Lagrange_1_2D();

  /// Destructor
  virtual ~dof_map_Lagrange_1_2D();

  /// Return a string identifying the dof map
  virtual const char* signature() const;

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(unsigned int d) const;

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m);

  /// Initialize dof map for given cell
  virtual void init_cell(const ufc::mesh& m,
                         const ufc::cell& c);

  /// Finish initialization of dof map for cells
  virtual void init_cell_finalize();

  /// Return the dimension of the global finite element function space
  virtual unsigned int global_dimension() const;

  /// Return the dimension of the local finite element function space
  virtual unsigned int local_dimension() const;

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const;

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   unsigned int facet) const;

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double** coordinates,
                                    const ufc::cell& c) const;

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const;

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const;

};


} // namespace

#endif
