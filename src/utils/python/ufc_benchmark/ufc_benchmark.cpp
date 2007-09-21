#include "ufc_benchmark.h"

int benchmark(const ufc::form & form, unsigned int n)
{
    ufc::mesh m;
    ufc::cell c;
    // FIXME: fill mesh and cell
    for(unsigned int i=0; i<n; i++)
    {
        //form->tabulate_tensor(...); // FIXME
    }
}
