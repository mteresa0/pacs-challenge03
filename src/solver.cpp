#include "solver.hpp"
#include "boundaries.hpp"
#include <iostream>
#include <cmath>
#include "omp.h"

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

    void Solver::evaluate_boundaries(std::vector<double> & u_local) const
    {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // update global lower
        if(rank==0)
            for (index_type i = 1; i<global_N-1;++i) {
                u_local[i] = bds.B_lower.fun(domain.a + i*domain.h);
            };
        
        
        index_type local_rows = (global_N%size > static_cast<unsigned int>(rank)) ?  global_N/size+1 : global_N/size;

        // update global upper
        if (rank==size-1){
            for (index_type i = 1; i<global_N-1;++i) {
                u_local[i+(local_rows-1)*global_N] = bds.B_upper.fun(domain.a + i*domain.h);
            }
        }

        for (index_type i = 0; i<local_rows; ++i)
        {
            // update right
            u_local[(i+1)*global_N-1] = bds.B_right.fun(get_global_row(rank, size, i));
            // update left
            u_local[i*global_N] = bds.B_left.fun(get_global_row(rank, size, i));
        }

    }

    std::vector<double> Solver::compute_solution() const
    {
        std::vector<double> u(domain.get_size_grid(), 0);

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int num_threads = omp_get_num_threads();

        const unsigned int min_local_size = global_N/size;
        const unsigned int remainder = global_N%size;

        index_type local_rows = (remainder > static_cast<unsigned int>(rank)) ?  min_local_size+1 : min_local_size;

        std::vector<int> start_index, local_size;
        start_index.resize(size); local_size.resize(size);
        int count = 0;
        for (unsigned int r = 0; r < static_cast<unsigned int>(size); ++r){
            local_size[r] = (remainder>r) ?  (min_local_size+1)*global_N : min_local_size*global_N;
            start_index[r] = count;
            count += local_size[r];
        }
        
        index_type global_starting_row = start_index[rank]/global_N;

        // std::cout << "Rank " << rank << " has " << local_rows << " rows\n";

        std::vector<double> old_local_u(global_N*local_rows, 0);

        evaluate_boundaries(old_local_u);

        std::vector<double> upper(global_N-2, 0);
        std::vector<double> lower(global_N-2, 0);

        std::vector<double> upper_to_send(global_N-2, 0);
        std::vector<double> lower_to_send(global_N-2, 0);

        bool converged = false;
        std::vector<double> new_local_u(old_local_u);

        std::vector<double> local_f(old_local_u);
        for (index_type j = 0; j<local_rows; ++j)
        {
            for (index_type i = 0; i<global_N; ++i)
                local_f[i+j*global_N] = f(domain.get_coord(i), domain.get_coord(j+global_starting_row));
        }

        double prod = 0;
        double err = 0;
        unsigned n;
        for(n = 0; n<max_it && !converged; ++n)
        {
            err = 0;
            if (local_rows>2){
                #pragma omp parallel for shared(new_local_u) num_threads(num_threads)
                for (index_type j = 1; j<(local_rows-1); ++j) 
                {
                    for (index_type i = 1; i<(global_N-1); ++i)
                    {
                        new_local_u[get_vector_index(i,j)] = (old_local_u[get_vector_index(i+1,j)] +
                                                                old_local_u[get_vector_index(i-1, j)] +
                                                                old_local_u[get_vector_index(i, j+1)] + 
                                                                old_local_u[get_vector_index(i, j-1)] + 
                                                                domain.h*domain.h*local_f[i+j*global_N])*0.25;
                        prod = new_local_u[get_vector_index(i,j)]-old_local_u[get_vector_index(i,j)];
                        err += prod*prod;
                    }
                }
            }

            if (rank!=0) {
                for (index_type i = 1; i<global_N-1; ++i)
                {
                    // j = 0;
                    new_local_u[get_vector_index(i,0)] = 0.25*(old_local_u[get_vector_index(i+1,0)] +
                                                            old_local_u[get_vector_index(i-1, 0)] +
                                                            old_local_u[get_vector_index(i, 1)] +
                                                            lower[i-1] +
                                                            domain.h*domain.h*local_f[i]);
                    prod = new_local_u[i]-old_local_u[i];
                    err += prod*prod;
                }
            }

            if (rank!=size-1) {
                for (index_type i = 1; i<global_N-1; ++i)
                {
                    index_type j = local_rows-1;
                    new_local_u[get_vector_index(i,j)] = 0.25*(old_local_u[get_vector_index(i+1,j)] +
                                                            old_local_u[get_vector_index(i-1, j)] +
                                                            old_local_u[get_vector_index(i, j-1)] +
                                                            upper[i-1] +
                                                            domain.h*domain.h*local_f[i+j*global_N]);
                    prod = new_local_u[get_vector_index(i,j)]-old_local_u[get_vector_index(i,j)];
                    err += prod*prod;
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            err *= domain.h;
            err = std::sqrt(err);

            converged = err<tol; 

            if (!converged)
            {
                MPI_Barrier(MPI_COMM_WORLD);

                MPI_Status status;
                // receive lower from rank-1
                if (rank!=0)
                {
                    MPI_Recv(lower.data(), global_N-2, MPI_DOUBLE, rank-1, (rank-1), MPI_COMM_WORLD, &status);
                }
                // sent local upper (for the upper rank is the lower)
                if (rank!=size-1)
                {
                    for (index_type i = 0; i<global_N-2; ++i){
                        upper_to_send[i] = new_local_u[1+i+(local_rows-1)*global_N];
                    }
                    MPI_Send(upper_to_send.data(), global_N-2, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
                }

                // sent local lower (to the lower rank is the upper)
                if (rank!=0)
                {
                    for (index_type i = 0; i<global_N-2; ++i){
                        lower_to_send[i] = new_local_u[i+1];
                    }
                    MPI_Send(lower_to_send.data(), global_N-2, MPI_DOUBLE, rank-1, (rank)+global_N, MPI_COMM_WORLD);
                }
                // receive upper data from rank+1
                if (rank!=size-1)
                {
                    MPI_Recv(upper.data(), global_N-2, MPI_DOUBLE, rank+1, (rank+1)+global_N, MPI_COMM_WORLD, &status);
                }

                old_local_u = new_local_u;
            }
            
        } // end for loop for computing

        if (rank == 0) std::cout << "conveged in " << n << " iterations\n"; 

        MPI_Allgatherv(new_local_u.data(), local_rows*global_N, MPI_DOUBLE, 
        u.data(), local_size.data(), start_index.data(), MPI_DOUBLE, MPI_COMM_WORLD);

        return u;
    }

} //  end namespace laplacian_solver