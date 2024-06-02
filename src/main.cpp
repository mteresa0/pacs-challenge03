#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#pragma GCC diagnostic pop

#include <iostream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <numbers>
#include "domain.hpp"
#include "solver.hpp"
#include "writeVTKfile.hpp"



int main(int argc, char ** argv) {
    
    using namespace laplacian_solver;

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (MPI_THREAD_MULTIPLE>provided)
        std::cerr << "MPI implementation do not provide enough support for MPI_THREAD_MULPTIPLE";

    constexpr double pi = M_PI;

    double a = 0.; double b = 1.;
    unsigned int N = std::stoul(argv[1]);

    Domain domain(N, a, b);

    source_type f = [](const double & x, const double & y) {return 8*pi*pi*sin(2*pi*x)*sin(2*pi*y);};
    source_type u_ex = [](const double & x, const double & y) {return sin(2*pi*x)*sin(2*pi*y);};

    Solver s1(domain, f, u_ex);

    std::vector<double> u_ = s1.compute_solution();

    if (rank == 0) write_VTK(domain, u_, "solution.vtk", "solution u");
    
    MPI_Finalize();

    return 0;

}