%module ClebschGordan

%{
#define SWIG_FILE_WITH_INIT
#include "ClebschGordan.h"
%}

namespace clebsch {
    const double EPS = 1e-12;

    // binomial coefficients
    class binomial_t {
        std::vector<int> cache;
        int N;

    public:
        int operator()(int n, int k);
    };
}
