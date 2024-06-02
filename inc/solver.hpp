#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include "domain.hpp"
#include <functional>
#include <mpi.h>
#include "boundaries.hpp"

namespace laplacian_solver {

    using index_type = std::size_t;
    using source_type=std::function<double(double, double)>;

    class Solver {
        const Domain domain;
        const Boundaries bds;
        const double tol;
        const unsigned max_it;

        const source_type f;
        const source_type u_ex;

    public:
        Solver(const Domain & _domain, const source_type & _f, const source_type & _u_ex, const double & _tol = 1e-6, const unsigned & _max_it = 1000) : 
        domain(_domain), f(_f), u_ex(_u_ex), tol(_tol), max_it(_max_it) ,bds(_domain) {};

        void evaluate_boundaries(std::vector<double> & );

        std::vector<double> compute_solution();

        inline void print() const { std::cout << domain.N << "\n";};

        inline index_type get_global_row(const int & rank, const int & size, const index_type & local_row) {
        return (domain.N%size>rank) ? rank*(1+domain.N/size) + local_row : (domain.N%size)+(domain.N/size)*rank + local_row;
        };     

        inline double f_discretized(index_type i, index_type j) {
            return f(domain.get_coord(i), domain.get_coord(j));
        };

        inline index_type get_vector_index(const index_type & i, const index_type & j) {
            return i+j*domain.N;
        }
    };


} // namespace laplacian_solver

#endif // SOLVER_HPP   