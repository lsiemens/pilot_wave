#ifndef MCINT_H
#define MCINT_H

#include <random>
#include <functional>
#include <iostream>

#include "testing/error_interval.h"

std::mt19937 get_rng() {
    std::random_device rd;
    std::seed_seq sequence{rd(), rd(), rd(), rd()};
    return std::mt19937(sequence);
}

ErrorInterval mc_integrate_uniform_1D(double a, double b, std::size_t N, std::function<double(double)> integrand) {
    auto rng = get_rng();
    auto uniform = std::uniform_real_distribution<double>(a, b);

    // welford's algorithm
    double mean_old = 0.;
    double mean = 0.;
    double M_2 = 0.;
    double v_i = 0.;
    for (std::size_t i = 1; i <= N; i++) {
        v_i = integrand(uniform(rng));
        mean_old = mean;
        mean += (v_i - mean)/static_cast<double>(i);
        M_2 += (v_i - mean_old)*(v_i - mean);
    }
    double variance = M_2/static_cast<double>(N);

    double value = (b - a)*mean;
    double error = (b - a)*std::sqrt(variance/static_cast<double>(N));
    return ErrorInterval(value, error);
}

double _riemann_sum_1D(double a, double b, std::size_t resolution, std::function<double(double)> integrand) {
    double dx = (b - a)/static_cast<double>(resolution);
    double value = 0;
    for (std::size_t i = 0; i < resolution; i++) {
        double x = a + dx*(static_cast<double>(i) + 0.5);
        value += dx*integrand(x);
    }
    return value;
}

ErrorInterval integrate_uniform_1D(double a, double b, std::size_t resolution, std::function<double(double)> integrand) {
    double I_half_resolution = _riemann_sum_1D(a, b, resolution/2, integrand);
    double I_full_resolution = _riemann_sum_1D(a, b, resolution, integrand);
    double error = 2*std::abs(I_half_resolution - I_full_resolution);
    return ErrorInterval(I_full_resolution, error);
}

double _riemann_sum_3D(double a, double b, std::size_t resolution, std::function<double(glm::dvec3)> integrand) {
    double dx = (b - a)/static_cast<double>(resolution);
    double dv = dx*dx*dx;
    double value = 0;
    for (std::size_t i = 0; i < resolution; i++) {
        double x = a + dx*(static_cast<double>(i) + 0.5);
        for (std::size_t j = 0; j < resolution; j++) {
            double y = a + dx*(static_cast<double>(j) + 0.5);
            for (std::size_t k = 0; k < resolution; k++) {
                double z = a + dx*(static_cast<double>(k) + 0.5);

                value += dv*integrand(glm::dvec3(x, y, z));
            }
        }
    }
    return value;
}

ErrorInterval integrate_uniform_3D(double a, double b, std::size_t resolution, std::function<double(glm::dvec3)> integrand) {
    double I_half_resolution = _riemann_sum_3D(a, b, resolution/2, integrand);
    double I_full_resolution = _riemann_sum_3D(a, b, resolution, integrand);
    double error = 2*std::abs(I_half_resolution - I_full_resolution);
    return ErrorInterval(I_full_resolution, error);
}

#endif
