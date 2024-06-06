#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <iostream>
#include <vector>

namespace laplacian_solver {
    using index_type = std::size_t;

    /// @brief describes a 2D square domain [a,b]x[a,b]
    /// each dimension is divided by N nodes
    struct Domain 
    {
        /// @brief number of nodes
        const index_type N;
        /// @brief bounds of the domain
        const double a, b;
        /// @brief distretization step
        const double h;

        // explicit Domain() : N(11), a(0.), b(1.), h(0.1) {};
        Domain(const index_type & _N, const double & _a, const double & _b) : 
        N(_N), a(_a), b(_b), h((_b-_a)/(_N-1)) {};

        /// @brief print domain data
        inline void print() const {
            std::cout << "Square domain [" << a << "," << b <<"]x[" << a << "," << b <<"]\n";
            std::cout << "N = " << N <<"\n";
            std::cout << "h = " << h <<"\n";
            return;
        };

        /// @brief get coordinates in either x or y since the two dimensions 
        /// have the same discretization and values
        /// @param i index (from 0 to N-1)
        /// @return i-th node in x or y direction
        inline double get_coord(const index_type & i) const {return h*i+a;};
        
    };

}; // namespace laplacian_solver

#endif // DOMAIN_HPP