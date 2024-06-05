#ifndef SOLVER_HPP
#define SOLVER_HPP

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#pragma GCC diagnostic pop

#include <vector>
#include "domain.hpp"
#include <functional>
#include "boundaries.hpp"

namespace laplacian_solver {

    using index_type = std::size_t;
    using source_type=std::function<double(double, double)>;

    class Solver {
        const Domain domain;
        const index_type global_N;
        const Boundaries bds;

        const source_type f;
        const source_type u_ex;

        const double tol;
        const unsigned max_it;

    public:
        Solver(const Domain & _domain, const source_type & _f, 
                const source_type & _u_ex, const double & _tol = 1e-6, 
                const unsigned & _max_it = 50000) : 
        domain(_domain), global_N(_domain.N), bds(_domain), f(_f), u_ex(_u_ex), tol(_tol), max_it(_max_it) {};

        void evaluate_boundaries(std::vector<double> & ) const;

        std::vector<double> compute_solution() const;

        double L2_norm(const std::vector<double> & ) const;

        inline void print() const { std::cout << global_N << "\n";};

        inline index_type get_global_row(const int & rank, const int & size, const index_type & local_row) const 
        {
        return (global_N%size > static_cast<unsigned int>(rank)) ? 
        rank*(1+global_N/size) + local_row : (global_N%size)+(global_N/size)*rank + local_row;
        };     

        inline double f_discretized(index_type i, index_type j) const {
            return f(domain.get_coord(i), domain.get_coord(j));
        };

        inline index_type get_vector_index(const index_type & i, const index_type & j) const {
            return (i+j*global_N);
        }

    };

} // namespace laplacian_solver

#endif // SOLVER_HPP   