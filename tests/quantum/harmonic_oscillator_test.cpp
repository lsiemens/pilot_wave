#include <iostream>

#include "core/math_util.h"
#include "quantum/harmonic_oscillator.h"
#include "testing/testing.h"

bool test_constants() {
    HarmonicOscillator ho(1.0);
    ho.set_energy_level(1);

    if (not is_close(ho.m_hbar, 1., "HarmonicOscillator.m_hbar")) {
        return false;
    }

    if (not is_close(ho.m_m_e, 1., "HarmonicOscillator.m_m_e")) {
        return false;
    }
  
    return true;
}

bool test_factorial(std::size_t N) {
    HarmonicOscillator ho(1.0);
    ho.set_energy_level((N + 1)*(N + 2)*(N + 3)/6);

    for (std::size_t n = 0; n <= N; n++) {
        double fac = 1.;
        for (std::size_t i = 1; i <= n; i++) {
            fac *= static_cast<double>(i);
        }

        std::string name = "HarmonicOscillator.factorial(" + std::to_string(n) + ")";
        if (not is_close(ho.factorial(n), fac, name)) {
            return false;
        }
    }

    return true;
}

double H_5(double x) {
    return 32.*std::pow(x, 5) - 160.*std::pow(x, 3) + 120.*x;
}

double H_10(double x) {
    return    1024.*std::pow(x, 10) - 23040.*std::pow(x, 8)
           +161280.*std::pow(x, 6) - 403200.*std::pow(x, 4)
           +302400.*std::pow(x, 2) - 30240.;
}

/// Test HarmonicOscillator.hermite_n against the 5th and 10th Hermite polynomial.
bool test_hermite() {
    HarmonicOscillator ho(1.0);

    double width = 5.;
    std::size_t resolution = 100;
    double dx = width/static_cast<double>(resolution);
    for (std::size_t i = 0; i <= resolution; i++) {
        double x = dx*static_cast<double>(i)-0.5*width;

        std::string name = "HarmonicOscillator.Hermite(5, " + std::to_string(x) + ")";
        if (not is_close(ho.hermite_n(5, x), H_5(x), name)) {
            return false;
        }

        name = "HarmonicOscillator.Hermite(10, " + std::to_string(x) + ")";
        if (not is_close(ho.hermite_n(10, x), H_10(x), name)) {
            return false;
        }
    }

    return true;
}

double psi_0_1D(double x, double omega) {
    double norm = std::pow(omega/PI_D, 1/4.);
    return norm*std::exp(-0.5*omega*x*x);
}

double psi_4_1D(double x, double omega) {
    double norm = std::pow(omega/PI_D, 1/4.)/std::sqrt(16.*24.);
    double H_4 = 16.*omega*omega*std::pow(x, 4) - 48.*omega*x*x + 12.;
    return norm*std::exp(-0.5*omega*x*x)*H_4;
}

/// Test HarmonicOscillator.psi_n_1D against the 0th and 4th state.
bool test_psi_1D(double omega) {
    HarmonicOscillator ho(omega);
    std::size_t N_max = 4;
    ho.set_energy_level((N_max + 1)*(N_max + 2)*(N_max + 3)/6);

    double width = 6/std::sqrt(2*omega);
    std::size_t resolution = 100;
    double dx = width/static_cast<double>(resolution);
    for (std::size_t i = 0; i <= resolution; i++) {
        double x = dx*static_cast<double>(i)-0.5*width;

        std::string name = "HarmonicOscillator.psi_1D(0, " + std::to_string(x) + ")";
        if (not is_close(ho.psi_n_1D(0, x), psi_0_1D(x, omega), name, 1e-5, 1e-4)) {
            return false;
        }

        name = "HarmonicOscillator.psi_1D(4, " + std::to_string(x) + ")";
        if (not is_close(ho.psi_n_1D(4, x), psi_4_1D(x, omega), name, 1e-5, 1e-4)) {
            return false;
        }
    }

    return true;
}

int main() {
    if (not test_constants()) {
        std::cerr << "\ttest_constants() failed." << std::endl;
        return 1;
    }

    if (not test_factorial(10)) {
        std::cerr << "\ttest_factorial(10) failed." << std::endl;
        return 1;
    }

    if (not test_hermite()) {
        std::cerr << "\ttest_hermite() failed." << std::endl;
        return 1;
    }

    if (not test_psi_1D(1.0)) {
        std::cerr << "\ttest_psi_1D(1.0) failed." << std::endl;
        return 1;
    }

    if (not test_psi_1D(1/3.)) {
        std::cerr << "\ttest_psi_1D(1/3) failed." << std::endl;
        return 1;
    }

    return 0;
}
