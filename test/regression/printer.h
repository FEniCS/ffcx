#ifndef PRINTER_H_INCLUDED
#define PRINTER_H_INCLUDED

#include <cassert>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>

class Printer
{
protected:

  // Precision in output of floats
  const std::size_t precision;
  const double epsilon;

  // Output stream
  std::ostream & os;

  // Indentation level
  int level;


public:

  Printer(std::ostream & os):
    precision(16),
    epsilon(1e-16),
    os(os),
    level(0)
  {}


  /// Indent to current level
  void indent()
  {
    for (int i=0; i<level; ++i)
      os << "    ";
  }

  /// Format name with optional integers appended
  std::string format_name(std::string name, int i=-1, int j=-1)
  {
    std::stringstream s;
    s << name;
    if (i >= 0)
      s << "_" << i;
    if (j >= 0)
      s << "_" << j;
    return s.str();
  }

  /// Format '<indent>"foo": ' properly
  void begin_entry(std::string name, int i=-1, int j=-1)
  {
    name = format_name(name, i, j);
    indent();
    os << '"' << name << '"' << ": ";
  }

  /// Begin an unnamed block
  void begin()
  {
    os << "{" << std::endl;
    ++level;
  }

  /// Begin a named block entry
  void begin(std::string name, int i=-1, int j=-1)
  {
    begin_entry(name, i, j);
    begin();
  }

  /// End current named or unnamed block
  void end()
  {
    --level;
    assert(level >= 0);
    indent();
    if (level > 0)
      os << "}," << std::endl;
    else
      os << "}" << std::endl;
  }

  /// Format a value type properly
  template<typename T>
  void print_value(T value);

  /// Set a named single value entry
  template<typename T>
  void print_scalar(std::string name, T value, int i=-1, int j=-1)
  {
    begin_entry(name, i, j);
    print_value(value);
    os << ", " << std::endl;
  }

  /// Set a named array valued entry
  template<typename T>
  void print_array(std::string name, int n, T * values, int i=-1, int j=-1)
  {
    begin_entry(name, i, j);
    os << "[";
    if (n > 0)
      print_value(values[0]);
    for (int k=1; k<n; ++k)
      {
        os << ", ";
        print_value(values[k]);
      }
    os << "], " << std::endl;
  }

  /// Set a named vector valued entry
  template<typename T>
  void print_vector(std::string name, typename std::vector<T> values, int i=-1, int j=-1)
  {
    begin_entry(name, i, j);
    os << "[";
    typename std::vector<T>::iterator k=values.begin();
    if (k!=values.end())
      {
        print_value(*k);
        ++k;
      }
    while (k!=values.end())
      {
        os << ", ";
        print_value(*k);
        ++k;
      }
    os << "], " << std::endl;
  }

};

/// Fallback formatting for any value type
template<typename T>
void Printer::print_value(T value)
{
  os << value;
}

/// Use precision for floats
template<>
void Printer::print_value(double value)
{
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(precision);
  if (std::abs(static_cast<double>(value)) < epsilon)
    os << "0.0";
  else
    os << value;
}

/// Use precision for floats
template<>
void Printer::print_value(float value)
{
  print_value(static_cast<double>(value));
}

/// Wrap strings in quotes
template<>
void Printer::print_value(std::string value)
{
  os << '"' << value << '"';
}

/// Wrap strings in quotes
template<>
void Printer::print_value(const char * value)
{
  os << '"' << value << '"';
}

#endif
