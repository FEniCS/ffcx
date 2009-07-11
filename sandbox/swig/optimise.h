#ifndef __OPTIMISE_H
#define __OPTIMISE_H

#import <string>
#import <vector>

class Symbol
{
public:

  Symbol(){}
  Symbol(std::string name, std::string type);

  Symbol copy();

  std::string name, type;
};

class Product
{
public:

  Product(std::vector<Symbol> vars);
//  Product(std::vector<int> vars);

  std::vector<Symbol> vars;
//  std::vector<int> vars;
};

#endif
