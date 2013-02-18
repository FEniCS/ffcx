
#include "simplejson.h"

#include <vector>
#include <map>
#include <fstream>

int main(int argc, char* argv[])
{
    std::ofstream af("a.output");
    JsonPrinter a(af);
    a.json_begin();
    a.json_begin("foo");
    a.json_print_scalar("a", 1);
    a.json_print_scalar("b", 2.3);
    std::vector<double> c1(2);
    c1[0] = 5.6;
    c1[1] = 6.5;
    a.json_print_vector("c", c1);
    a.json_begin("bar");
    a.json_print_scalar("a", 2);
    a.json_print_scalar("b", 3.4);
    a.json_end();
    a.json_begin("only_a");
    a.json_print_scalar("avaluei", 7);
    a.json_print_scalar("avaluef", "3.14159", 1);
    a.json_print_scalar("avalues", "somename", 1, 2);
    double aarr[3] = { 1.2, 3.4, 5.6 };
    a.json_print_array("avaluearr", 3, aarr);
    a.json_end();
    a.json_end();
    a.json_end();

    std::ofstream bf("b.output");
    JsonPrinter b(bf);
    b.json_begin();
    b.json_begin("foo");
    b.json_print_scalar("a", 1);
    b.json_print_scalar("b", 2.3);
    std::vector<double> c2;
    c2.push_back(5.6);
    c2.push_back(6.51);
    b.json_print_vector("c", c2);
    b.json_begin("bar");
    b.json_print_scalar("a", 2);
    b.json_print_scalar("b", 5.4);
    b.json_end();
    b.json_begin("only_b");
    b.json_print_scalar("bvaluei", 7);
    b.json_print_scalar("bvaluef", 3.14159);
    std::string someothername("someothername");
    b.json_print_scalar("bvalues", someothername);
    double barr[3] = { 1.2, 3.4, 5.6 };
    b.json_print_array("bvaluearr", 3, barr, 4, 5);
    b.json_end();
    b.json_end();
    b.json_end();

    return 0;
}
