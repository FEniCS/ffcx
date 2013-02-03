
#include <vector>
#include <list>
#include <map>

#include <cassert>
#include <string>
#include <iostream>
#include <fstream>

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

    /// Set a named single variable entry
    template<typename T>
    void operator()(const std::string & name, T value)
    {
        begin_entry(name);
        output_formatted(value);
        out << ", " << std::endl;
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


int main(int argc, char* argv[])
{
    std::ofstream af("a.output");
    SimpleJsonModel a(af);
    a.begin();
    a.begin("foo");
    a("a", 1);
    a("b", 2.3);
    a.begin("bar");
    a("a", 2);
    a("b", 3.4);
    a.end();
    a.begin("only_a");
    a("avaluei", 7);
    a("avaluef", "3.14159");
    a("avalues", "somename");
    a.end();
    a.end();
    a.end();
 
    std::ofstream bf("b.output");
    SimpleJsonModel b(bf);
    b.begin();
    b.begin("foo");
    b("a", 1);
    b("b", 2.3);
    b.begin("bar");
    b("a", 2);
    b("b", 5.4);
    b.end();
    b.begin("only_b");
    b("bvaluei", 7);
    b("bvaluef", 3.14159);
    std::string someothername("someothername");
    b("bvalues", someothername);
    b.end();
    b.end();
    b.end();

    return 0;
}
