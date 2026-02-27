#include "core/math_util.h"

#include <vector>
#include <cmath>
#include <iostream>

static std::vector<double> cumsum_ln = {0.0};

inline void update_factorial(std::size_t n) {
    std::size_t n_old = cumsum_ln.size() - 1;
    if (n > n_old) {
        cumsum_ln.resize(n + 1);
        for (std::size_t i = n_old + 1; i < n + 1; i++) {
            cumsum_ln[i] = std::log(i) + cumsum_ln[i -  1];
        }
    }
}

double ln_factorial(std::size_t n) {
    update_factorial(n);
    return cumsum_ln[n];
}

double factorial(std::size_t n) {
    update_factorial(n);
    return std::exp(cumsum_ln[n]);
}
