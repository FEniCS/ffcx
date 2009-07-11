%module optimise

%{
#include <vector>
%}

%{
#include "optimise.h"
%}

%include std_string.i
%include std_vector.i

%template(vector_int)     std::vector<int>;
%typedef std::vector<int> vector_int;

%template(vector_symbol)     std::vector<Symbol>;
%typedef std::vector<Symbol> vector_symbol;

%include "optimise.h"


