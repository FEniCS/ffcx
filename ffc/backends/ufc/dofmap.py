# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2017.

ufc_dofmap_declaration = """
extern "C" ufc_dofmap* create_{factory_name}();
"""


ufc_dofmap_combined = """
struct {factory_name}
{{
};
"""


dofmap_combined = """
class %(classname)s: public ufc::dofmap
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::dofmap()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const char * signature() const final override
  {
%(signature)s
  }

  int64_t num_global_support_dofs() const final override
  {
%(num_global_support_dofs)s
  }

  int64_t num_element_support_dofs() const final override
  {
%(num_element_support_dofs)s
  }

  int64_t num_element_dofs() const final override
  {
%(num_element_dofs)s
  }

  int64_t num_facet_dofs() const final override
  {
%(num_facet_dofs)s
  }

  int64_t num_entity_dofs(int64_t d) const final override
  {
%(num_entity_dofs)s
  }

  int64_t num_entity_closure_dofs(int64_t d) const final override
  {
%(num_entity_closure_dofs)s
  }

  void tabulate_dofs(int64_t * dofs,
                     const int64_t* num_global_entities,
                     const int64_t** entity_indices) const final override
  {
%(tabulate_dofs)s
  }

  void tabulate_facet_dofs(int64_t * dofs,
                           int64_t facet) const final override
  {
%(tabulate_facet_dofs)s
  }

  void tabulate_entity_dofs(int64_t * dofs,
                            int64_t d, int64_t i) const final override
  {
%(tabulate_entity_dofs)s
  }

  void tabulate_entity_closure_dofs(int64_t * dofs,
                                    int64_t d, int64_t i) const final override
  {
%(tabulate_entity_closure_dofs)s
  }


  int64_t num_sub_dofmaps() const final override
  {
%(num_sub_dofmaps)s
  }

  ufc::dofmap * create_sub_dofmap(int64_t i) const final override
  {
%(create_sub_dofmap)s
  }

  ufc::dofmap * create() const final override
  {
%(create)s
  }

};
"""

dofmap_header = """
class %(classname)s: public ufc::dofmap
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const char * signature() const final override;

  int64_t num_global_support_dofs() const final override;

  int64_t num_element_support_dofs() const final override;

  int64_t num_element_dofs() const final override;

  int64_t num_facet_dofs() const final override;

  int64_t num_entity_dofs(int64_t d) const final override;

  int64_t num_entity_closure_dofs(int64_t d) const final override;

  void tabulate_dofs(int64_t * dofs,
                     const int64_t* num_global_entities,
                     const int64_t** entity_indices) const final override;

  void tabulate_facet_dofs(int64_t * dofs,
                           int64_t facet) const final override;

  void tabulate_entity_dofs(int64_t * dofs,
                            int64_t d, int64_t i) const final override;

  void tabulate_entity_closure_dofs(int64_t * dofs,
                            int64_t d, int64_t i) const final override;

  int64_t num_sub_dofmaps() const final override;

  ufc::dofmap * create_sub_dofmap(int64_t i) const final override;

  ufc::dofmap * create() const final override;

};
"""

dofmap_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::dofmap()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const char * %(classname)s::signature() const
{
%(signature)s
}

int64_t %(classname)s::num_global_support_dofs() const
{
%(num_global_support_dofs)s
}

int64_t %(classname)s::num_element_support_dofs() const
{
%(num_element_support_dofs)s
}

int64_t %(classname)s::num_element_dofs() const
{
%(num_element_dofs)s
}

int64_t %(classname)s::num_facet_dofs() const
{
%(num_facet_dofs)s
}

int64_t %(classname)s::num_entity_dofs(int64_t d) const
{
%(num_entity_dofs)s
}

int64_t %(classname)s::num_entity_closure_dofs(int64_t d) const
{
%(num_entity_closure_dofs)s
}

void %(classname)s::tabulate_dofs(int64_t * dofs,
                                  const int64_t* num_global_entities,
                                  const int64_t** entity_indices) const
{
%(tabulate_dofs)s
}

void %(classname)s::tabulate_facet_dofs(int64_t * dofs,
                                        int64_t facet) const
{
%(tabulate_facet_dofs)s
}

void %(classname)s::tabulate_entity_dofs(int64_t * dofs,
                                         int64_t d, int64_t i) const
{
%(tabulate_entity_dofs)s
}

void %(classname)s::tabulate_entity_closure_dofs(int64_t * dofs,
                                             int64_t d, int64_t i) const
{
%(tabulate_entity_closure_dofs)s
}

int64_t %(classname)s::num_sub_dofmaps() const
{
%(num_sub_dofmaps)s
}

ufc::dofmap * %(classname)s::create_sub_dofmap(int64_t i) const
{
%(create_sub_dofmap)s
}

ufc::dofmap * %(classname)s::create() const
{
%(create)s
}
"""
