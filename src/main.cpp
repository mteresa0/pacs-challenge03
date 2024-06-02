
#include <iostream>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <numbers>
#include "domain.hpp"
#include "solver.hpp"

int main(int argc, char ** argv) {
    // using laplacian_solver;

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    constexpr double pi = M_PI;

    laplacian_solver::Domain domain(11, 0.,1.);

    laplacian_solver::source_type f = [](const double & x, const double & y) {return 8*pi*pi*sin(2*pi*x)*sin(2*pi*y);};
    laplacian_solver::source_type u_ex = [](const double & x, const double & y) {return sin(2*pi*x)*sin(2*pi*y);};

    std::vector<double> u(domain.get_size_grid(),0);

    laplacian_solver::Solver s1(domain, f, u_ex);
    s1.print();

    MPI_Finalize();

    // auto f = [](const double & x, const double & y) {return 8*pi*pi*sin(2*pi*x)*sin(2*pi*y);};
    // auto u_ex = [](const double & x, const double & y) {return sin(2*pi*x)*sin(2*pi*y);};

    // double a_x, b_x, a_y, b_y;
    // a_x = 0; b_x = 1;
    // a_y = 0; b_y = 1;
    // index_type N_x = 11;
    // index_type N_y = 11;
    // double h_x = (b_x-a_x)/N_x;
    // double h_y = (b_y-a_y)/N_y;

    // std::vector<double> U(N_x*N_y,0);



    // if (MPI_THREAD_MULTIPLE>provided)
    //     std::cerr << "MPI implementation do not provide enough support for MPI_THREAD_MULPTIPLE";

    // int rank,size;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);

    // #pragma omp parallel 
    // {
    // std::string str = "hi from thread " + std::to_string(omp_get_thread_num()) + " of process " + std::to_string(rank) + "\n";
    // std::cout << str;
    // };

    // MPI_Finalize();

    return 0;

}