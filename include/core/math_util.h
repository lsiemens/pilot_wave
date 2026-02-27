#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include <numbers>

constexpr double PI_D = std::numbers::pi_v<double>;
constexpr float PI_F = std::numbers::pi_v<float>;

double ln_factorial(std::size_t n);

double factorial(std::size_t n);

#endif
