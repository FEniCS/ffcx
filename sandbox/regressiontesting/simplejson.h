#ifndef SIMPLE_JSON_H__INCLUDED
#define SIMPLE_JSON_H__INCLUDED

#include <cassert>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>

/// A simplified json data model, including only hierarchial dicts of single values
class SimpleJsonModel
{
public:
    SimpleJsonModel(std::ostream & out):
        level(0),
        out(out)
    {
    }

    /// Begin an unnamed block
    void begin()
    {
        out << "{" << std::endl;
        ++level;
    }

    /// Begin a named block entry
    void begin(const std::string & name)
    {
        begin_entry(name);
        begin();
    }

    /// Set a named single value entry
    template<typename T>
    void print_scalar(std::string name, T value, int i=-1, int j=-1)
    {
        std::stringstream s;
        s << name;
        if (i >= 0) s << "_" << i;
        if (j >= 0) s << "_" << j;
        name = s.str();

        begin_entry(name);
        print_formatted(value);
        out << ", " << std::endl;
    }

    /// Set a named array valued entry
    template<typename T>
    void print_array(std::string name, int n, T * values, int i=-1, int j=-1)
    {
        std::stringstream s;
        s << name;
        if (i >= 0) s << "_" << i;
        if (j >= 0) s << "_" << j;
        name = s.str();

        begin_entry(name);
        out << "[";
        if (n > 0)
          print_formatted(values[0]);
        for (int k=1; k<n; ++k)
        {
            out << ", ";
            print_formatted(values[k]);
        }
        out << "], " << std::endl;
    }

    /// Set a named vector valued entry
    template<typename T>
    void print_vector(std::string name, typename std::vector<T> values, int i=-1, int j=-1)
    {
        std::stringstream s;
        s << name;
        if (i >= 0) s << "_" << i;
        if (j >= 0) s << "_" << j;
        name = s.str();

        begin_entry(name);
        out << "[";
        typename std::vector<T>::iterator k=values.begin();
        if (k!=values.end())
        {
            print_formatted(*k);
            ++k;
        }
        while (k!=values.end())
        {
            out << ", ";
            print_formatted(*k);
            ++k;
        }
        out << "], " << std::endl;
    }

    /// Set a named list valued entry
    template<typename T>
    void print_list(std::string name, typename std::list<T> values, int i=-1, int j=-1)
    {
        std::stringstream s;
        s << name;
        if (i >= 0) s << "_" << i;
        if (j >= 0) s << "_" << j;
        name = s.str();

        begin_entry(name);
        out << "[";
        typename std::list<T>::iterator k=values.begin();
        if (k!=values.end())
        {
            print_formatted(*k);
            ++k;
        }
        while (k!=values.end())
        {
            out << ", ";
            print_formatted(*k);
            ++k;
        }
        out << "], " << std::endl;
    }

    /// End current named or unnamed block
    void end()
    {
        --level;
        assert(level >= 0);
        indent();
        if (level > 0)
            out << "}," << std::endl;
        else
            out << "}" << std::endl;
    }

protected:

    int level;
    std::ostream & out;

    /// Indent to current level
    void indent()
    {
        for (int i=0; i<level; ++i)
            out << "    ";
    }

    /// Output line with indention and newline
    void output(const std::string & line)
    {
        indent();
        out << line << std::endl;
    }

    /// Format a value type properly
    template<typename T>
    void print_formatted(T value);

    /// Format '<indent>"foo": ' properly
    void begin_entry(const std::string & name)
    {
        indent();
        out << '"' << name << '"' << ": ";
    }
};

/// Fallback formatting for any value type
template<typename T>
void SimpleJsonModel::print_formatted(T value)
{
    out << value;
}

/// Wrap strings in quotes
template<>
void SimpleJsonModel::print_formatted(std::string value)
{
    out << '"' << value << '"';
}

/// Wrap strings in quotes
template<>
void SimpleJsonModel::print_formatted(const char * value)
{
    out << '"' << value << '"';
}

#endif
