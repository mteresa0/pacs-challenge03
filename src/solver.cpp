#include "solver.hpp"
#include "boundaries.hpp"
#include <iostream>
#include <cmath>

namespace laplacian_solver{
    void Boundaries::check_boundaries_compatibility(const double & a, const double & b)
    {
        // lower - right
        if (B_lower.fun(b) != B_right.fun(a))
            std::cerr << "boundary conditions are not compatibile (lower - right)";

                
        // upper - right
        if (B_upper.fun(b) != B_right.fun(b))
            std::cerr << "boundary conditions are not compatibile (upper - right)";

                
        // lower - left
        if (B_lower.fun(a) != B_left.fun(a))
            std::cerr << "boundary conditions are not compatibile (lower - left)";

                
        // upper - right
        if (B_upper.fun(a) != B_left.fun(b))
            std::cerr << "boundary conditions are not compatibile (upper - left)";

        return;
    };

    void Solver::evaluate_boundaries(std::vector<double> & u_local)
    {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // update global upper
        if(rank==0)
            for (index_type i = 1; i<domain.N-1;++i) {
                u_local[i] = bds.B_lower.fun(domain.a + i*domain.h);
            };
        
        
        index_type local_rows = (domain.N%size>rank) ?  domain.N/size+1 : domain.N/size;

        // update global upper
        if (rank==size-1){
            for (index_type i = 1; i<domain.N-1;++i) {
                u_local[i+(local_rows-1)*domain.N] = bds.B_upper.fun(domain.a + i*domain.h);
            }
        }

        for (index_type i = 0; i<local_rows; ++i)
        {
            // update right
            u_local[(i+1)*domain.N-1] = bds.B_right.fun(get_global_row(rank, size, i));
            // update left
            u_local[i*domain.N] = bds.B_left.fun(get_global_row(rank, size, i));
        }

    }

    std::vector<double> Solver::compute_solution()
    {
        std::vector<double> u(domain.get_size_grid(), 0);

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        index_type local_rows = (domain.N%size>rank) ?  domain.N/size+1 : domain.N/size;

        std::vector<double> old_local_u(domain.N*local_rows, 0);

        evaluate_boundaries(old_local_u);

        std::vector<double> upper(domain.N-2, 0);
        std::vector<double> lower(domain.N-2, 0);

        bool converged = false;
        std::vector<double> new_local_u(old_local_u);

        double prod = 0;
        double err = 0;
        for(unsigned n = 0; n<max_it && !converged; ++n)
        {
            err = 0;
            for (index_type i = 1; i<domain.N-1; ++i){
                for (index_type j = 1; j>local_rows-1; ++j)
                {
                    new_local_u[get_vector_index(i,j)] = 1/4*(old_local_u[get_vector_index(i+1,j)] +
                                                            old_local_u[get_vector_index(i-1, j)] +
                                                            old_local_u[get_vector_index(i, j+1)] +
                                                            old_local_u[get_vector_index(i, j-1)] +
                                                            domain.h*domain.h*f_discretized(i,get_global_row(rank, size, j)));
                    prod = new_local_u[get_vector_index(i,j)]-old_local_u[get_vector_index(i,j)];
                    err += prod*prod;
                }
            }

            if (rank!=0) {
                for (index_type i = 1; i<domain.N-1; ++i)
                {
                    new_local_u[get_vector_index(i,0)] = 1/4*(old_local_u[get_vector_index(i+1,0)] +
                                                            old_local_u[get_vector_index(i-1, 0)] +
                                                            old_local_u[get_vector_index(i, 1)] +
                                                            lower[i-1] +
                                                            domain.h*domain.h*f_discretized(i,get_global_row(rank, size, 0)));
                    prod = new_local_u[i]-old_local_u[i];
                    err += prod*prod;
                }
            }

            if (rank!=size-1) {
                for (index_type i = 1; i<domain.N-1; ++i)
                {
                    index_type j = local_rows-1;
                    new_local_u[get_vector_index(i,j)] = 1/4*(old_local_u[get_vector_index(i+1,j)] +
                                                            old_local_u[get_vector_index(i-1, j)] +
                                                            old_local_u[get_vector_index(i, j-1)] +
                                                            upper[i-1] +
                                                            domain.h*domain.h*f_discretized(i,get_global_row(rank, size, 0)));
                    prod = new_local_u[get_vector_index(i,j)]-old_local_u[get_vector_index(i,j)];
                    err += prod*prod;
                }
            }

            err *= domain.h;
            err = std::sqrt(err);

            converged = err<tol;

            if (!converged)
            {
                
            }
            
        }


        return u;
    }

} //  end namespace laplacian_solver