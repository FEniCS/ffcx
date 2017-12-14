# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2017

function_combined = """
class %(classname)s: public ufc::function
{%(members)s
public:

  %(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::function()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  void evaluate(double * values,
                const double * coordinates,
                const ufc::cell& c) const final override
  {
%(evaluate)s
  }

};
"""

function_header = """
class %(classname)s: public ufc::function
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  void evaluate(double * values,
                const double * coordinates,
                const ufc::cell& c) const final override;

};
"""

function_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::function()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

void %(classname)s::evaluate(double * values,
                             const double * coordinates,
                             const ufc::cell& c) const
{
%(evaluate)s
}
"""
