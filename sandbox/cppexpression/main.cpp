
#include <iostream>
#include <dolfin.h>

namespace uflacs
{
    namespace test_uflfiles_test
    {
        class Expression_g: public dolfin::Expression
        {
        public:
            virtual void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
            {
                // Implementation of: ((x)[0]) * ((x)[1])
                double s[1];
                s[0] = x[0] * x[1]; // ((x)[0]) * ((x)[1])
                values[0] = s[0];
            }
        };
    }
}

int main()
{
  uflacs::test::Expression_g g;
  return 0;
}
