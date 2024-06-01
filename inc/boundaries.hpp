
#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP

#include <functional>
#include <vector>

namespace laplacian_solver
{
        
    enum boundary_edge {
        upper,
        right,
        lower,
        left
    };

    template <boundary_edge BE>
    class Boundary
    {
        const unsigned int N;
    public:
        Boundary(const unsigned int & _N) : N(_N) {};
        virtual ~Boundary() = default;

        inline unsigned int get_size() {return N;};
        std::vector<unsigned int> get_boundary_indexes();

        virtual void evaluate_boundary(std::vector<double> & u) const {return;};
    };

    template<boundary_edge BE>
    class Dirichlet_Boundary : public Boundary {
        std::function(double(double)) fun;
    public:
        Dirichlet_Boundary(const unsigned int & _N, std::function(double(double)) _fun = [](const double &){return 0;}) : Boundary(_N), fun(_fun) {};
        void evaluate_boundary(std::vector<double> & u) const override;
    };

    template<boundary_edge BE>
    std::vector<unsigned int> Boundary<BE>::get_boundary_indexes()
    {
        return std::vector<unsigned int>();
    }

    // template <typename T>
    // std::vector<unsigned int> Boundary<T, upper>::get_boundary_indexes() {
        
    //     std::vector<unsigned int> idx(N,0);

    //     for (unsigned int i = 0; i<N; ++i)
    //     {
    //         idx[i] = N*(N-1) + i;
    //     }

    //     return idx;
    // }


} // namespace laplacian_solver

#endif // BOUNDARIES_HPP

