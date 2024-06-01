#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <iostream>

namespace laplacian_solver {
    class Domain 
    {
        const unsigned int N;
        const double a, b;
        const double h;
    public:
        Domain() = default;
        Domain(const unsigned int & _N, const double & _a, const double & _b) : 
        N(_N), a(_a), b(_b), h((_b-_a)/(_N-1)) {};

        inline void print() const {
            std::cout << "Square domain [" << a << "," << b <<"]x[" << a << "," << b <<"]\n";
            std::cout << "N = " << N <<"\n";
            std::cout << "h = " << h <<"\n";
        };

        friend class Solver;
    };

}; // namespace laplacian_solver

#endif // DOMAIN_HPP