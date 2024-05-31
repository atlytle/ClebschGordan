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

    // Eq. (19) and (25)
    class weight {
        std::vector<int> elem;

    public:
        // the N in "SU(N)"
        const int N;

        // create a non-initialized weight
        weight(int N);

        // create irrep weight of given index
        // Eq. (C2)
        weight(int N, int index);

        // assign from another instance
        clebsch::weight &operator=(const clebsch::weight &w);

        // access elements of this weight (k = 1, ..., N)
        int &operator()(int k);
        const int &operator()(int k) const;

        // allow python access to data through swig
        %rename(__getitem__) operator[];
        %rename(__setitem__) operator[];
        %extend {
            int __getitem__(int k) {
                return (*($self))[k];
            }

            void __setitem__(int k, int v) {
                (*($self))[k] = v;
            }
        }


        // compare weights
        // Eq. (C1)
        bool operator<(const weight &w) const;
        bool operator==(const weight &w) const;

        // element-wise sum of weights
        clebsch::weight operator+(const weight &w) const;

        // returns the index of this irrep weight (index = 0, 1, ...)
        // Eq. (C2)
        int index() const;

        // returns the dimension of this irrep weight
        // Eq. (22)
        long long dimension() const;
    };
class coefficients {
        std::map<std::vector<int>, double> clzx;

        // access Clebsch-Gordan coefficients in convenient manner
        void set(int factor1_state,
                 int factor2_state,
                 int multiplicity_index,
                 int irrep_state,
                 double value);

        // internal functions, doing most of the work
        void highest_weight_normal_form(); // Eq. (37)
        void compute_highest_weight_coeffs(); // Eq. (36)
        void compute_lower_weight_coeffs(int multip_index, int state, std::vector<char> &done); // Eq. (40)

    public:
        // the N in "SU(N)"
        const int N;

        // save irreps and their dimensions for later use
        const weight factor1, factor2, irrep;
        const int factor1_dimension, factor2_dimension, irrep_dimension;

        // outer multiplicity of irrep in this decomposition
        const int multiplicity;

        // construct all Clebsch-Gordan coefficients of this decomposition
        coefficients(const weight &irrep, const weight &factor1, const weight &factor2);

        // access Clebsch-Gordan coefficients (read-only)
        // multiplicity_index = 0, ..., multiplicity - 1
        // factor1_state = 0, ..., factor1_dimension - 1
        // factor2_state = 0, ..., factor2_dimension - 1
        // irrep_state = 0, ..., irrep_dimension
        double operator()(int factor1_state,
                          int factor2_state,
                          int multiplicity_index,
                          int irrep_state) const;
    };
}
