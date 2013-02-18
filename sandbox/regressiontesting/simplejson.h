#ifndef SIMPLE_JSON_H__INCLUDED
#define SIMPLE_JSON_H__INCLUDED

#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>

/// A simplified json data model, including only hierarchial dicts of single values
class JsonPrinter
{
protected:
  // Precision in output of floats
  const std::size_t precision;
  const double epsilon;

protected:
  std::ostream & json_os;
  int level;

  /// Indent to current level
  void json_indent()
  {
    for (int i=0; i<level; ++i)
      json_os << "    ";
  }

  /// Format a value type properly
  template<typename T>
  void json_print_value(T value);

  /// Format '<indent>"foo": ' properly
  void json_begin_entry(const std::string & name)
  {
    json_indent();
    json_os << '"' << name << '"' << ": ";
  }

public:
  JsonPrinter(std::ostream & json_os):
    precision(16),
    epsilon(1e-16),

    json_os(json_os),
    level(0)
  {
  }

  /// Begin an unnamed block
  void json_begin()
  {
    json_os << "{" << std::endl;
    ++level;
  }

  /// Begin a named block entry
  void json_begin(const std::string & name)
  {
    json_begin_entry(name);
    json_begin();
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

  /// Set a named single value entry
  template<typename T>
  void json_print_scalar(std::string name, T value, int i=-1, int j=-1)
  {
    name = json_format_name(name, i, j);
    json_begin_entry(name);
    json_print_value(value);
    json_os << ", " << std::endl;
  }

  /// Set a named array valued entry
  template<typename T>
  void json_print_array(std::string name, int n, T * values, int i=-1, int j=-1)
  {
    name = json_format_name(name, i, j);
    json_begin_entry(name);
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
    name = json_format_name(name, i, j);
    json_begin_entry(name);
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
};

/// Fallback formatting for any value type
template<typename T>
void JsonPrinter::json_print_value(T value)
{
  json_os << value;
}

/// Use precision for floats
template<>
void JsonPrinter::json_print_value(double value)
{
  json_os.precision(precision);
  if (std::abs(static_cast<double>(value)) < epsilon)
    json_os << "0.0";
  else
    json_os << value;
}

/// Use precision for floats
template<>
void JsonPrinter::json_print_value(float value)
{
  json_os.precision(precision);
  if (std::abs(static_cast<double>(value)) < epsilon)
    json_os << "0.0";
  else
    json_os << value;
}

/// Wrap strings in quotes
template<>
void JsonPrinter::json_print_value(std::string value)
{
  json_os << '"' << value << '"';
}

/// Wrap strings in quotes
template<>
void JsonPrinter::json_print_value(const char * value)
{
  json_os << '"' << value << '"';
}

#endif
