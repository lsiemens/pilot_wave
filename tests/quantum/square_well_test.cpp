#include <iostream>

#include "quantum/square_well.h"
#include "testing/testing.h"

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

bool test_wave_function(std::size_t n_max, std::size_t resolution=10) {
    double width = 1.0;

    SquareWell sw(width);
    sw.set_energy_level(n_max);

    for (std::size_t n = 0; n < n_max; n++) {
        double dx = width/static_cast<double>(resolution);
        double dv = dx*dx*dx;
        double integral = 0;
        for (std::size_t i = 0; i < resolution; i++) {
            double x = dx*(static_cast<double>(i) + 0.5);
            for (std::size_t j = 0; j < resolution; j++) {
                double y = dx*(static_cast<double>(j) + 0.5);
                for (std::size_t k = 0; k < resolution; k++) {
                    double z = dx*(static_cast<double>(k) + 0.5);
                    integral += dv*std::norm(sw.psi_n(glm::dvec3(x, y, z), n));
                }
            }
        }

        std::string name = "Integral of |psi_" + std::to_string(n) + "(x)|^2";
        if (not is_close(integral, 1., name)) {
            return false;
        } 
    }

    for (std::size_t n = 1; n < n_max; n++) {
        for (std::size_t m = 0; m <= n; m++) {
                
            double dx = width/static_cast<double>(resolution);
            double dv = dx*dx*dx;
            double integral = 0;
            for (std::size_t i1 = 0; i1 < resolution; i1++) {
                double x1 = dx*(static_cast<double>(i1) + 0.5);
                for (std::size_t j1 = 0; j1 < resolution; j1++) {
                    double y1 = dx*(static_cast<double>(j1) + 0.5);
                    for (std::size_t k1 = 0; k1 < resolution; k1++) {
                        double z1 = dx*(static_cast<double>(k1) + 0.5);
                        glm::dvec3 x_vec(x1, y1, z1);

                        integral += dv*sw.psi_n(x_vec, n)*sw.psi_n(x_vec, m);
                    }
                }
            }

            std::string name = "Integral of |psi_" + std::to_string(n) + "(x)"
                                           "*psi_" + std::to_string(m) + "(x)|^2";
            if (n == m) {
                if (not is_close(integral, 1., name)) {
                    return false;
                } 
            } else {
                if (not is_close(integral, 0., name)) {
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

    if (not test_wave_function(5)) {
        std::cerr << "\ttest_wave_function(5) failed." << std::endl;
        return 1;
    }

    return 0;
}
