#include "optimise.h"
//#include <iostream>

Symbol::Symbol(std::string name, std::string type) :
name(name), type(type)
{
//  std::cout << "Initialise Symbol" << std::endl;
}

Symbol Symbol::copy()
{
//  return Symbol(name, count, type);
  return *this;
}

Product::Product(std::vector<Symbol> vars) : vars(vars)
//Product::Product(std::vector<int> vars) : vars(vars)
{
}


