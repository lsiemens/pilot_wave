#include <iostream>

#include "quantum/square_well.h"

#include "testing/testing.h"
#include "testing/integrate.h"

bool test_constants() {
    SquareWell sw(1.0);
    sw.set_energy_level(1);

    if (not is_close(sw.m_hbar, 1., "SquareWell.m_hbar")) {
        return false;
    }

    if (not is_close(sw.m_m_e, 1., "SquareWell.m_m_e")) {
        return false;
    }
  
    return true;
}

bool test_orthonorm_3D(double omega) {
    SquareWell sq(omega);
    std::size_t N_max = 10;
    sq.set_energy_level(N_max);
    double x_max = 10.;

    std::size_t resolution = 100;
    for (std::size_t n = 0; n < N_max; n++) {
        for (std::size_t m = 0; m <= n; m++) {

            auto integrand_3D = [&sq, n, m](glm::dvec3 x) {
                return sq.psi_n(x, n)*sq.psi_n(x, m);
            };

            auto value = integrate_uniform_3D(-x_max, x_max, resolution, integrand_3D);

            std::string name = "integral(-10, 10) |psi_3D(" + std::to_string(n) + ", x)psi_3D(" + std::to_string(m) + ", x)|^2";
            if (n == m) {
                if (not is_close(value, 1., name)) {
                    return false;
                }
            } else {
                if (not is_close(value, 0., name)) {
                    return false;
                }
            }
        }
    }
    return true;
}

int main() {
    if (not test_constants()) {
        std::cerr << "\ttest_constants() failed." << std::endl;
        return 1;
    }

    if (not test_orthonorm_3D(1.0)) {
        std::cerr << "\ttest_orthonorm_3D(1.0) failed." << std::endl;
        return 1;
    }

    return 0;
}
