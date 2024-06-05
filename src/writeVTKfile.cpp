#include "writeVTKfile.hpp"

namespace laplacian_solver::write_output{

    void write_VTK_2D_domain(const unsigned int &Nx, const unsigned int Ny, const double &hx, const double &hy, const std::vector<double> &data, const std::string &filename, const std::string &title)
    {
        std::ofstream vtk_file(filename);

        if (!vtk_file.is_open()) {
            std::cerr << "error in opening VTK file \n";
        }

        vtk_file << "# vtk DataFile Version 3.0\n";
        vtk_file << title <<"\nASCII\n";

        vtk_file << "DATASET STRUCTURED_POINTS\n";
        vtk_file << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
        vtk_file << "ORIGIN 0.0 0.0 0.0\n";
        vtk_file << "SPACING " << hx << " " << hy << " 1\n";
        if (data.size()!=Nx*Ny) std::cerr << "data size do not match dimensions of the grid";
        vtk_file << "POINT_DATA " << data.size() <<"\n";
        vtk_file << "SCALARS scalars double\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (unsigned int i = 0; i<data.size(); ++i){
            vtk_file << data[i] << std::endl;
        };

        vtk_file.close();

        return;
    }

    void write_convergence(const std::vector<unsigned int> &nodes, const std::vector<double> &spacings, const std::vector<double> &norms, const std::string &filename)
    {
        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "error in opening CSV file " << filename << " \n";
        };

        unsigned int size = nodes.size();

        file << "N, h, L2_norm" << std::endl;

        for (unsigned int i = 0; i<size; ++i){
            file << nodes[i] << ", "<< spacings[i]<< ", " << norms[i] << std::endl;
        }

        file.close();

        return;
    };
}