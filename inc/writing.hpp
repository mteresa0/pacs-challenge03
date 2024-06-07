#ifndef WRITE_VTK_2D_HPP
#define WRITE_VTK_2D_HPP

#include <iostream>
#include <fstream>
#include "domain.hpp"

namespace laplacian_solver::write_output{

    void write_VTK_2D_domain(const unsigned int & Nx, const unsigned int Ny, const double & hx, const double & hy,
                const std::vector<double> & data, const std::string & filename, const std::string & title = "data");

    inline void write_VTK_2D_domain(const unsigned int & N, const double & h, const std::vector<double> & data, 
    const std::string & filename, const std::string & title = "data")
    {
        write_VTK_2D_domain(N, N, h, h, data, filename, title);
        return;
    }

    inline void write_VTK_2D_domain(const laplacian_solver::Domain & domain ,const std::vector<double> & data, 
    const std::string & filename, const std::string & title = "data")
    {
        write_VTK_2D_domain(domain.N, domain.h, data, filename, title);
        return;
    };

    void write_convergence(const std::vector<unsigned int> & nodes, const std::vector<double> & spacings, const std::vector<double> & norms , const std::string & filename);

} // end namespace 

#endif // WRITE_VTK_2D_HPP
