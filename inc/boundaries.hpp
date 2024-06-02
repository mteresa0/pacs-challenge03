
#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP

#include <functional>
#include <vector>

namespace laplacian_solver
{
    using index_type = std::size_t;
    using boundary_function = std::function<double(const double &)>;
        
    enum boundary_edge {
        upper,
        right,
        lower,
        left
    };

    template <boundary_edge BE>
    struct Boundary
    {
        const index_type N;

        boundary_function fun;

        Boundary(const index_type & _N, const boundary_function & _fun = [](const double &){return 0;}) : 
        N(_N), fun(_fun) {};

        /// @brief 
        /// @param  
        /// @return 
        index_type get_boundary_index(const index_type &) const;
    };

    struct Boundaries
    {
        const Boundary<upper> B_upper;
        const Boundary<lower> B_lower;
        const Boundary<left> B_left;
        const Boundary<right> B_right;

        void check_boundaries_compatibility(const double & a, const double & b);

        Boundaries(const Domain & domain): B_upper(Boundary<upper>(domain.N)), B_lower(Boundary<lower>(domain.N)), 
        B_left(Boundary<left>(domain.N)), B_right(Boundary<right>(domain.N)) {check_boundaries_compatibility(domain.a, domain.b);};

        Boundaries(const Boundary<upper> & _upper, const Boundary<lower> & _lower, 
        const Boundary<left> & _left, const Boundary<right> & _right) : 
        B_upper(_upper), B_lower(_lower), B_left(_left), B_right(_right) {};
    };

    template<>
    inline index_type Boundary<upper>::get_boundary_index(const index_type & k) const
    {return (N-1)*N + k;};

    template<>
    inline index_type Boundary<lower>::get_boundary_index(const index_type & k) const
    {return k;};

    template<>
    inline index_type Boundary<left>::get_boundary_index(const index_type & k) const
    {return k*(N);};

    template<>
    inline index_type Boundary<right>::get_boundary_index(const index_type & k) const
    {return (k+1)*N - 1;};

} // namespace laplacian_solver

#endif // BOUNDARIES_HPP

