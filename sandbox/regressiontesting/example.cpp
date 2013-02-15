
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
    a.output_scalar("a", 1);
    a.output_scalar("b", 2.3);
    std::vector<double> c1(2);
    c1[0] = 5.6;
    c1[1] = 6.5;
    a.output_array("c", c1);
    a.begin("bar");
    a.output_scalar("a", 2);
    a.output_scalar("b", 3.4);
    a.end();
    a.begin("only_a");
    a.output_scalar("avaluei", 7);
    a.output_scalar("avaluef", "3.14159");
    a.output_scalar("avalues", "somename");
    a.end();
    a.end();
    a.end();
 
    std::ofstream bf("b.output");
    SimpleJsonModel b(bf);
    b.begin();
    b.begin("foo");
    b.output_scalar("a", 1);
    b.output_scalar("b", 2.3);
    std::list<double> c2;
    c2.push_back(5.6);
    c2.push_back(6.51);
    b.output_array("c", c2);
    b.begin("bar");
    b.output_scalar("a", 2);
    b.output_scalar("b", 5.4);
    b.end();
    b.begin("only_b");
    b.output_scalar("bvaluei", 7);
    b.output_scalar("bvaluef", 3.14159);
    std::string someothername("someothername");
    b.output_scalar("bvalues", someothername);
    b.end();
    b.end();
    b.end();

    return 0;
}
