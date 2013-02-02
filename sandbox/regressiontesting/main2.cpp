
#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>

class DataModel
{
public:
    DataModel(std::ostream & out):
        level(0),
        out(out)
    {
    }

    void indent()
    {
        for (int i=0; i<level; ++i)
            out << "    ";
    }

    void output(const std::string & line)
    {
        indent();
        out << line << std::endl;
    }

    void begin_entry(const std::string & name)
    {
        indent();
        out << '"' << name << '"' << ": ";
    }

    void begin_outer_block()
    {
        out << "{" << std::endl;
        ++level;
    }

    void begin_block(const std::string & name)
    {
        begin_entry(name);
        begin_outer_block();
    }

    template<typename T>
    void entry(const std::string & name, const T & value)
    {
        begin_entry(name);
        out << value << ", " << std::endl;
    }

    void end_block()
    {
        --level;
        indent();
        out << "}," << std::endl;
    }

    void end_outer_block()
    {
        --level;
        indent();
        out << "}" << std::endl;
    }

    int level;
    std::ostream & out;
};

int main()
{
    DataModel model(std::cout);

    model.output("#!/usr/bin/env python");
    model.output("");

    model.output("a = (");
    model.begin_outer_block();
    model.begin_block("foo");
    model.entry("a", 1);
    model.entry("b", 2.3);
    model.begin_block("bar");
    model.entry("a", 2);
    model.entry("b", 3.4);
    model.end_block();
    model.end_block();
    model.end_outer_block();
    model.output(")");
    model.output("");

    model.output("b = (");
    model.begin_outer_block();
    model.begin_block("foo");
    model.entry("a", 1);
    model.entry("b", 2.3);
    model.begin_block("bar");
    model.entry("a", 2);
    model.entry("b", 5.4);
    model.end_block();
    model.end_block();
    model.end_outer_block();
    model.output(")");
    model.output("");

    model.output("from recdiff import *");
    model.output("print_recdiff(recdiff(a,b))");
    model.output("");

    return 0;
}
