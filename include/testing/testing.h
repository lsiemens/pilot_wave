#ifndef TESTING_H
#define TESTING_H

#include <string>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <iomanip>

#include "testing/error_interval.h"

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
        std::cerr << "\tabs(ε) = [ " << diff << " ]" << std::endl;
        return false;
    }

    return true;
}

inline bool is_close(ErrorInterval a, double b, std::string name="", double min_error=1e-8, double max_rerror=1e-5) {

    double rerror = a.m_error/std::abs(std::max(a.m_value, min_error));
    if ((rerror > max_rerror) and (a.m_value > min_error)) {
        std::cerr << std::setprecision(5) << "Failed: " << name << " [ " << a.m_value << " ± " << std::max(a.m_error, min_error) << " ] != " << b << " is_close failed!" << std::endl;
        std::cerr << "\trelative error = [ " << rerror << " ]" << std::endl;
        return false;
    }
    double diff  = std::abs(a.m_value - b);
    double error = std::max(a.m_error, min_error);
    if (diff > error) {
        int sigfig = static_cast<int>(1 - std::floor(std::log10(error)));

        std::cerr << std::setprecision(sigfig) << "Failed: " << name << " [ " << a.m_value << " ± " << error << " ] != " << b << " is_close failed!" << std::endl;
        std::cerr << "\trange = [ " << a.m_value - error << ", " << a.m_value + error << " ]" << std::endl;
        return false;
    }

    return true;
}

#endif
