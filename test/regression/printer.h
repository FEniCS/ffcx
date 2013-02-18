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


protected:
  ////// Counter based printing code:

  // Output stream
  std::ostream & os;
  // Global counter for results
  std::size_t counter;

  // Format name with counter and optional integers appended
  std::string counted_format_name(std::string name, int i=-1, int j=-1)
  {
    std::stringstream s;
    s << counter++ << "_";
    s << name;
    if (i >= 0) s << "_" << i;
    if (j >= 0) s << "_" << j;
    return s.str();
  }

  // Function for printing a single value
  template <class value_type>
  void counted_print_value(value_type value)
  {
    os.precision(precision);
    if (std::abs(static_cast<double>(value)) < epsilon)
      os << "0";
    else
      os << value;
  }

  // Function for beginning named result block
  void counted_begin(std::string name)
  {
    os << std::endl;
    os << "Testing " << name << std::endl;
    os << "----------------------" << std::endl;
  }

  // Function for printing scalar result
  template <class value_type>
  void counted_print_scalar(std::string name, value_type value, int i=-1, int j=-1)
  {
    name = counted_format_name(name, i, j);
    os << name << " = ";
    counted_print_value(value);
    os << std::endl;
  }

  // Function for printing array result
  template <class value_type>
  void counted_print_array(std::string name, unsigned int n, value_type* values, int i=-1, int j=-1)
  {
    name = counted_format_name(name, i, j);
    os << name << " =";
    for (std::size_t i = 0; i < n; i++)
      {
        os << " ";
        counted_print_value(values[i]);
      }
    os << std::endl;
  }

  // Function for printing vector result
  template <class value_type>
  void counted_print_vector(std::string name, std::vector<value_type> values, int i=-1, int j=-1)
  {
    name = counted_format_name(name, i, j);
    os << name << " =";
    for (std::size_t i = 0; i < values.size(); i++)
      {
        os << " ";
        counted_print_value(values[i]);
      }
    os << std::endl;
  }


protected:
  ////// Json based printing code:

  // Json output stream
  std::ostream & json_os;
  // Indentation level
  int level;


  /// Indent to current level
  void json_indent()
  {
    for (int i=0; i<level; ++i)
      json_os << "    ";
  }

  /// Format name with optional integers appended
  std::string json_format_name(std::string name, int i=-1, int j=-1)
  {
    std::stringstream s;
    s << name;
    if (i >= 0) s << "_" << i;
    if (j >= 0) s << "_" << j;
    return s.str();
  }

  /// Format a value type properly
  template<typename T>
  void json_print_value(T value);

  /// Format '<indent>"foo": ' properly
  void json_begin_entry(std::string name, int i=-1, int j=-1)
  {
    name = json_format_name(name, i, j);
    json_indent();
    json_os << '"' << name << '"' << ": ";
  }

  /// Begin an unnamed block
  void json_begin()
  {
    json_os << "{" << std::endl;
    ++level;
  }

  /// Begin a named block entry
  void json_begin(std::string name, int i=-1, int j=-1)
  {
    json_begin_entry(name, i, j);
    json_begin();
  }

  /// Set a named single value entry
  template<typename T>
  void json_print_scalar(std::string name, T value, int i=-1, int j=-1)
  {
    json_begin_entry(name, i, j);
    json_print_value(value);
    json_os << ", " << std::endl;
  }

  /// Set a named array valued entry
  template<typename T>
  void json_print_array(std::string name, int n, T * values, int i=-1, int j=-1)
  {
    json_begin_entry(name, i, j);
    json_os << "[";
    if (n > 0)
      json_print_value(values[0]);
    for (int k=1; k<n; ++k)
      {
        json_os << ", ";
        json_print_value(values[k]);
      }
    json_os << "], " << std::endl;
  }

  /// Set a named vector valued entry
  template<typename T>
  void json_print_vector(std::string name, typename std::vector<T> values, int i=-1, int j=-1)
  {
    json_begin_entry(name, i, j);
    json_os << "[";
    typename std::vector<T>::iterator k=values.begin();
    if (k!=values.end())
      {
        json_print_value(*k);
        ++k;
      }
    while (k!=values.end())
      {
        json_os << ", ";
        json_print_value(*k);
        ++k;
      }
    json_os << "], " << std::endl;
  }

  /// End current named or unnamed block
  void json_end()
  {
    --level;
    assert(level >= 0);
    json_indent();
    if (level > 0)
      json_os << "}," << std::endl;
    else
      json_os << "}" << std::endl;
  }

public:
  Printer(std::ostream & os, std::ostream & json_os):
    precision(16),
    epsilon(1e-16),

    os(os),
    counter(0),

    json_os(json_os),
    level(0)
  {}

  // Function for beginning unnamed result block
  void begin()
  {
    json_begin();
  }

  // Function for beginning named result block
  void begin(std::string name, int i=-1, int j=-1)
  {
    json_begin(name, i, j);
    counted_begin(name);
  }

  // Function for ending named result block
  void end()
  {
    json_end();
  }

  // Function for printing scalar result
  template <class value_type>
  void print_scalar(std::string name, value_type value, int i=-1, int j=-1)
  {
    json_print_scalar(name, value, i, j);
    counted_print_scalar(name, value, i, j);
  }

  // Function for printing array result
  template <class value_type>
  void print_array(std::string name, unsigned int n, value_type* values, int i=-1, int j=-1)
  {
    json_print_array(name, n, values, i, j);
    counted_print_array(name, n, values, i, j);
  }

  // Function for printing vector result
  template <class value_type>
  void print_vector(std::string name, std::vector<value_type> values, int i=-1, int j=-1)
  {
    json_print_vector(name, values, i, j);
    counted_print_vector(name, values, i, j);
  }
};

/// Fallback formatting for any value type
template<typename T>
void Printer::json_print_value(T value)
{
  json_os << value;
}

/// Use precision for floats
template<>
void Printer::json_print_value(double value)
{
  json_os.precision(precision);
  if (std::abs(static_cast<double>(value)) < epsilon)
    json_os << "0.0";
  else
    json_os << value;
}

/// Use precision for floats
template<>
void Printer::json_print_value(float value)
{
  json_os.precision(precision);
  if (std::abs(static_cast<double>(value)) < epsilon)
    json_os << "0.0";
  else
    json_os << value;
}

/// Wrap strings in quotes
template<>
void Printer::json_print_value(std::string value)
{
  json_os << '"' << value << '"';
}

/// Wrap strings in quotes
template<>
void Printer::json_print_value(const char * value)
{
  json_os << '"' << value << '"';
}

#endif
