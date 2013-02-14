
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
