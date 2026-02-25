#ifndef TESTING_H
#define TESTING_H

#include <string>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <iomanip>

#include "testing/error_interval.h"
#include <gtest/gtest.h>

inline bool is_close(const double& a, const double& b,
                     double rtol=1e-5, double atol=1e-8) {
    double diff  = std::abs(a - b);
    return diff <= (atol + rtol*diff);
}

inline bool is_close(const ErrorInterval& a, const double& b,
                     double min_error=1e-8, double max_rerror=1e-5) {
    double rerror = a.m_error/std::abs(std::max(a.m_value, min_error));
    if ((rerror > max_rerror) and (a.m_value > min_error)) {
        return false;
    }

    double diff  = std::abs(a.m_value - b);
    double error = std::max(a.m_error, min_error);
    return diff <= error;
}

inline ::testing::AssertionResult IsClose(const char* expr_a, const char* expr_b,
                                          const double& a, const double& b,
                                          double rtol=1e-5, double atol=1e-8) {
    if (is_close(a, b, rtol, atol)) {
        return ::testing::AssertionSuccess();
    }

    std::ostringstream oss;
    oss << std::setprecision(10)
        << expr_a << " and " << expr_b
        << " are not close: "
        << a << " != " << b
        << " (rtol=" << rtol << ", atol=" << atol << ")";
    return ::testing::AssertionFailure() << oss.str();
}

inline ::testing::AssertionResult IsClose(const char* expr_a, const char* expr_b,
                                          const ErrorInterval& a, const double& b,
                                          double min_error=1e-8, double max_rerror=1e-5) {
    if (is_close(a, b, min_error, max_rerror)) {
        return ::testing::AssertionSuccess();
    }

    double error = std::max(a.m_error, min_error);
    int sigfig = static_cast<int>(1 - std::floor(std::log10(error)));
    std::ostringstream oss;
    oss << std::setprecision(sigfig)
        << expr_a << " and " << expr_b
        << " are not close: "
        << "[ " << a.m_value << " ± " << a.m_error << "] != " << b
        << " (min_error=" << min_error << ", max_rerror=" << max_rerror << ")";
    return ::testing::AssertionFailure() << oss.str();
}

#endif
