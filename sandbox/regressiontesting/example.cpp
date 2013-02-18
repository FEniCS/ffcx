
#include "simplejson.h"

#include <vector>
#include <list>
#include <map>
#include <fstream>

int main(int argc, char* argv[])
{
    std::ofstream af("a.output");
    SimpleJsonModel a(af);
    a.begin();
    a.begin("foo");
    a.print_scalar("a", 1);
    a.print_scalar("b", 2.3);
    std::vector<double> c1(2);
    c1[0] = 5.6;
    c1[1] = 6.5;
    a.print_vector("c", c1);
    a.begin("bar");
    a.print_scalar("a", 2);
    a.print_scalar("b", 3.4);
    a.end();
    a.begin("only_a");
    a.print_scalar("avaluei", 7);
    a.print_scalar("avaluef", "3.14159", 1);
    a.print_scalar("avalues", "somename", 1, 2);
    double aarr[3] = { 1.2, 3.4, 5.6 };
    a.print_array("avaluearr", 3, aarr);
    a.end();
    a.end();
    a.end();
 
    std::ofstream bf("b.output");
    SimpleJsonModel b(bf);
    b.begin();
    b.begin("foo");
    b.print_scalar("a", 1);
    b.print_scalar("b", 2.3);
    std::list<double> c2;
    c2.push_back(5.6);
    c2.push_back(6.51);
    b.print_list("c", c2);
    b.begin("bar");
    b.print_scalar("a", 2);
    b.print_scalar("b", 5.4);
    b.end();
    b.begin("only_b");
    b.print_scalar("bvaluei", 7);
    b.print_scalar("bvaluef", 3.14159);
    std::string someothername("someothername");
    b.print_scalar("bvalues", someothername);
    double barr[3] = { 1.2, 3.4, 5.6 };
    b.print_array("bvaluearr", 3, barr, 4, 5);
    b.end();
    b.end();
    b.end();

    return 0;
}
