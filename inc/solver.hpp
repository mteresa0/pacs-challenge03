#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include "domain.hpp"
#include <functional>
// #include "boundaries.hpp"

namespace laplacian_solver {

    using source_type=std::function<double(double, double)>;

    class Solver {
        const Domain domain;
        // Boundary<upper> B_upper;
        // Boundary<lower> B_lower;
        // Boundary<left> B_left;
        // Boundary<right> B_right;

        const source_type f;
        const source_type u_ex;

        std::vector<double> u;

    public:
        Solver(const Domain & _domain, const source_type & _f, 
        const source_type & _u_ex, const std::vector<double> & _u)
         : domain(_domain), f(_f), u_ex(_u_ex), u(_u) {};


        inline void print() const { std::cout << domain.N;};
        
    };


} // namespace laplacian_solver

#endif // SOLVER_HPP   