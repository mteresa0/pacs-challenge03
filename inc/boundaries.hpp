
#ifndef BOUNDARIES_HPP
#define BOUNDARIES_HPP

#include <functional>
#include <vector>

namespace laplacian_solver
{
    using boundary_function = std::function<double(const double &)>;

    struct Boundaries
    {
        // boundary function for each side
        const boundary_function B_upper;
        const boundary_function B_right;
        const boundary_function B_lower;
        const boundary_function B_left;

        void check_boundaries_compatibility(const double & a, const double & b);

        Boundaries(const double & a, const double & b, const boundary_function & bup, 
        const boundary_function & bri, const boundary_function & blo, 
        const boundary_function & ble): 
        B_upper(bup), B_right(bri), B_lower(blo), B_left(ble)
        {check_boundaries_compatibility(a, b);};
        
        explicit Boundaries(): 
        B_upper([](const double &){return 0;}),
        B_right([](const double &){return 0;}), 
        B_lower([](const double &){return 0;}),
        B_left ([](const double &){return 0;}) {};
        

    };

} // namespace laplacian_solver

#endif // BOUNDARIES_HPP

