#ifndef WRITE_VTK_HPP
#define WRITE_VTK_HPP

#include <iostream>
#include <fstream>
#include "domain.hpp"

void write_VTK(const unsigned int & N, const double & h, const std::vector<double> & data, const std::string & filename, const std::string & title = "data")
{
    std::ofstream vtk_file(filename);

    if (!vtk_file.is_open()) {
        std::cerr << "error in opening VTK file \n";
    }

    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << title <<"\nASCII\n";

    vtk_file << "DATASET STRUCTURED_POINTS\n";
    vtk_file << "DIMENSIONS " << N << " " << N << " 1\n";
    vtk_file << "ORIGIN 0.0 0.0 0.0\n";
    vtk_file << "SPACING " << h << " " << h << " 1\n";
    vtk_file << "POINT_DATA " << data.size() <<"\n";
    vtk_file << "SCALARS scalars double\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (unsigned int i = 0; i<data.size(); ++i){
        vtk_file << data[i] << std::endl;
    };

    vtk_file.close();

    return;
};

void write_VTK(const laplacian_solver::Domain & domain ,const std::vector<double> & data, const std::string & filename, const std::string & title = "data"){
    write_VTK(domain.N, domain.h, data, filename, title);
};

#endif // WRITE_VTK_HPP