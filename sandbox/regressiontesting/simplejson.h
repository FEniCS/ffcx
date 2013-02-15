#ifndef SIMPLE_JSON_H__INCLUDED
#define SIMPLE_JSON_H__INCLUDED

#include <cassert>
#include <string>
#include <vector>
#include <list>
#include <iostream>

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
    void output_scalar(const std::string & name, T value)
    {
        begin_entry(name);
        output_formatted(value);
        out << ", " << std::endl;
    }

    /// Set a named vector valued entry
    template<typename T>
    void output_array(const std::string & name, typename std::vector<T> values)
    {
        begin_entry(name);
        out << "[";
        typename std::vector<T>::iterator i=values.begin();
        if (i!=values.end())
        {
            output_formatted(*i);
            ++i;
        }
        while (i!=values.end())
        {
            out << ", ";
            output_formatted(*i);
            ++i;
        }
        out << "], " << std::endl;
    }

    /// Set a named list valued entry
    template<typename T>
    void output_array(const std::string & name, typename std::list<T> values)
    {
        begin_entry(name);
        out << "[";
        typename std::list<T>::iterator i=values.begin();
        if (i!=values.end())
        {
            output_formatted(*i);
            ++i;
        }
        while (i!=values.end())
        {
            out << ", ";
            output_formatted(*i);
            ++i;
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
    void output_formatted(T value);

    /// Format '<indent>"foo": ' properly
    void begin_entry(const std::string & name)
    {
        indent();
        out << '"' << name << '"' << ": ";
    }
};

/// Fallback formatting for any value type
template<typename T>
void SimpleJsonModel::output_formatted(T value)
{
    out << value;
}

/// Wrap strings in quotes
template<>
void SimpleJsonModel::output_formatted(std::string value)
{
    out << '"' << value << '"';
}

/// Wrap strings in quotes
template<>
void SimpleJsonModel::output_formatted(const char * value)
{
    out << '"' << value << '"';
}

#endif
