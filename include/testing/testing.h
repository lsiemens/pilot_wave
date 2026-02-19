#ifndef TESTING_H
#define TESTING_H

#include <string>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <iomanip>

template <typename T, typename U>
inline bool is_close(T a, U b, std::string name="",
                     std::common_type_t<T, U> rtol=1e-5,
                     std::common_type_t<T, U> atol=1e-8) {
    using V = std::common_type_t<T, U>;

    static_assert(std::is_floating_point_v<V>, "Must be a floting point type!");

    V a_v = static_cast<V>(a);
    V b_v = static_cast<V>(b);

    V diff  = std::abs(a_v - b_v);
    if (not (diff <= (atol + rtol*diff))) {
        std::cerr << std::setprecision(10) << "Failed: " << name << " [ " << a << " ] != " << b << " is_close failed!" << std::endl;
        return false;
    }

    return true;
}

#endif
