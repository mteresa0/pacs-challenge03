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
    using source_type=std::function<double(double, double)>;

    class Solver {
        /// @brief struct containing the domain parameters
        const Domain domain;
        /// @brief number of nodes for each dimension
        const unsigned int global_N;
        /// @brief struct containing boundaries
        const Boundaries bds;

        /// @brief source function
        const source_type f;
        /// @brief exact solution
        const source_type u_ex;

        /// @brief tollerance step per rank
        const double tol;
        /// @brief max iterations
        const unsigned int max_it;

    public:
        Solver(const Domain & _domain, const source_type & _f, 
                const source_type & _u_ex, const double & _tol = 1e-6, 
                const unsigned & _max_it = 50000, Boundaries _bds = (Boundaries())) : 
        domain(_domain), global_N(_domain.N), bds(_bds), f(_f), u_ex(_u_ex), tol(_tol), max_it(_max_it) {};

        void evaluate_boundaries(std::vector<double> & ) const;

        std::vector<double> compute_solution(const std::string & ) const;

        double L2_norm(const std::vector<double> & ) const;

        /// @brief return vector index considering domain.N columns
        /// @param i 
        /// @param j
        inline unsigned int get_vector_index(const unsigned int & i, const unsigned int & j) const {
            return (i+j*global_N);
        }

    };

} // namespace laplacian_solver

#endif // SOLVER_HPP   