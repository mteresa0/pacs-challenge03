#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#include <omp.h>
#pragma GCC diagnostic pop

#include <iostream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <numbers>
#include "domain.hpp"
#include "solver.hpp"
#include "writing.hpp"

int main(int argc, char ** argv) {
    
    using namespace laplacian_solver;

    int provided;
    int req = MPI_THREAD_MULTIPLE;
    MPI_Init_thread(&argc, &argv, req, &provided);
    
    if (req>provided)
        std::cerr << "MPI implementation do not provide enough support for MPI_THREAD_MULPTIPLE";
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    constexpr double pi = M_PI;
    source_type f = [](const double & x, const double & y) {return 8*pi*pi*sin(2*pi*x)*sin(2*pi*y);};
    source_type u_ex = [](const double & x, const double & y) {return sin(2*pi*x)*sin(2*pi*y);};

    double a = 0.; double b = 1.; 

    Boundaries bds;

    unsigned int num_threads = 0;

    #pragma omp parallel
    {
    num_threads = omp_get_num_threads();
    }

    std::string convergence_filename = "outputs/convergence_"
                                    + std::to_string(size) 
                                    + "_" + std::to_string(num_threads) 
                                    +".csv";

    std::vector<unsigned int> nodes; 
    nodes.resize(argc-2);
    std::vector<double> norms_l2; 
    norms_l2.resize(argc-2);
    std::vector<double> spacings; 
    spacings.resize(argc-2);

    std::string times_filename = argv[1];

    for (unsigned int i = 0; i < static_cast<unsigned int>(argc)-2; ++i)
    {
        nodes[i] = std::stoul(argv[i+2]);
        Domain domain(nodes[i], a, b);

        Solver s1(domain, f, u_ex, 1e-9, 1e5, bds);

        std::vector<double> u_ = s1.compute_solution(times_filename);

        norms_l2[i] = (s1.L2_norm(u_));
        spacings[i] = (domain.h);

        if (rank == 0) {
            write_output::write_VTK_2D_domain(domain, u_, "outputs/solution"+std::to_string(nodes[i])+".vtk", "solution u");            
        }
        
    }

    if (rank == 0) write_output::write_convergence(nodes, spacings, norms_l2, convergence_filename);
    
    MPI_Finalize();

    return 0;

}