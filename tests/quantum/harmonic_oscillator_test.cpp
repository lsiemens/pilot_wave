#include <iostream>

#include "core/math_util.h"
#include "quantum/harmonic_oscillator.h"

#include "testing/testing.h"
#include "testing/integrate.h"
#include <gtest/gtest.h>

double H_5(double x) {
    return 32.*std::pow(x, 5) - 160.*std::pow(x, 3) + 120.*x;
}

double H_10(double x) {
    return    1024.*std::pow(x, 10) - 23040.*std::pow(x, 8)
           +161280.*std::pow(x, 6) - 403200.*std::pow(x, 4)
           +302400.*std::pow(x, 2) - 30240.;
}

/// Test HarmonicOscillator.hermite_n against the 5th and 10th Hermite polynomial.
TEST(HarmonicOscillatorTest, Hermite) {
    HarmonicOscillator ho(1.0);

    double width = 5.;
    std::size_t resolution = 100;
    double dx = width/static_cast<double>(resolution);
    for (std::size_t i = 0; i <= resolution; i++) {
        double x = dx*static_cast<double>(i)-0.5*width;

        std::string name = "HarmonicOscillator.Hermite(5, " + std::to_string(x) + ")";
        EXPECT_PRED_FORMAT2(IsClose, ho.hermite_n(5, x), H_5(x)) << name;

        name = "HarmonicOscillator.Hermite(10, " + std::to_string(x) + ")";
        EXPECT_PRED_FORMAT2(IsClose, ho.hermite_n(10, x), H_10(x)) << name;
    }
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
TEST(HarmonicOscillatorTest, Wavefunction1D) {
    double omega = 1.0;
    HarmonicOscillator ho(omega);
    std::size_t N_max = 4;
    ho.set_energy_level_index((N_max + 1)*(N_max + 2)*(N_max + 3)/6);

    double width = 6/std::sqrt(2*omega);
    std::size_t resolution = 100;
    double dx = width/static_cast<double>(resolution);
    for (std::size_t i = 0; i <= resolution; i++) {
        double x = dx*static_cast<double>(i)-0.5*width;

        std::string name = "HarmonicOscillator.psi_1D(0, " + std::to_string(x) + ")";
        EXPECT_PRED_FORMAT2(IsClose, ho.psi_n_1D(0, x), psi_0_1D(x, omega)) << name;

        name = "HarmonicOscillator.psi_1D(4, " + std::to_string(x) + ")";
        EXPECT_PRED_FORMAT2(IsClose, ho.psi_n_1D(4, x), psi_4_1D(x, omega)) << name;
    }
}

/// Test HarmonicOscillator.dpsi_n_1D_dx against the derivative of the 0th and 4th state.
TEST(HarmonicOscillatorTest, WavefunctionGradient1D) {
    double omega = 1.0;
    HarmonicOscillator ho(omega);
    std::size_t N_max = 4;
    ho.set_energy_level_index((N_max + 1)*(N_max + 2)*(N_max + 3)/6);

    double width = 6/std::sqrt(2*omega);
    std::size_t resolution = 100;
    double dx = width/static_cast<double>(resolution);
    double delta_x = 1e-5; // offset used in the central difference method
    for (std::size_t i = 0; i <= resolution; i++) {
        double x = dx*static_cast<double>(i)-0.5*width;

        std::string name = "HarmonicOscillator.dpsi_1D_dx(0, " + std::to_string(x) + ")";
        double dpsi_1D_dx = (psi_0_1D(x + delta_x, omega) - psi_0_1D(x - delta_x, omega))/(2*delta_x);
        EXPECT_PRED_FORMAT2(IsClose, ho.dpsi_n_1D_dx(0, x), dpsi_1D_dx) << name;

        name = "HarmonicOscillator.dpsi_1D_dx(4, " + std::to_string(x) + ")";
        dpsi_1D_dx = (psi_4_1D(x + delta_x, omega) - psi_4_1D(x - delta_x, omega))/(2*delta_x);
        EXPECT_PRED_FORMAT2(IsClose, ho.dpsi_n_1D_dx(4, x), dpsi_1D_dx) << name;
    }
}

TEST(HarmonicOscillatorTest, Orthonorm1D) {
    double omega = 1.0;
    HarmonicOscillator ho(omega);
    std::size_t N_max = 10;
    ho.set_energy_level_index((N_max + 1)*(N_max + 2)*(N_max + 3)/6);
    double x_max = 10.;

    std::size_t resolution = 100;
    for (std::size_t n = 0; n < N_max; n++) {
        for (std::size_t m = 0; m <= n; m++) {

            auto integrand_1D = [&ho, n, m](double x) {
                return ho.psi_n_1D(n, x)*ho.psi_n_1D(m, x);
            };

            auto value = integrate_uniform_1D(-x_max, x_max, resolution, integrand_1D);

            std::string name = "integral1D(-10, 10) |psi_1D(" + std::to_string(n) + ", x)psi_1D(" + std::to_string(m) + ", x)|^2";
            if (n == m) {
                EXPECT_PRED_FORMAT2(IsClose, value, 1.) << name;
            } else {
                EXPECT_PRED_FORMAT2(IsClose, value, 0.) << name;
            }
        }
    }
}

struct QN {
    int n_x, n_y, n_z;

    std::vector<int> to_vec() {
        return {n_x, n_y, n_z};
    }
    std::string to_str(std::size_t index) const {
        std::string expected_str = "Energy level [" + std::to_string(index) + "] quantum numbers: (" + std::to_string(n_x) + "," + std::to_string(n_y) + "," + std::to_string(n_z) + ")";
        return expected_str;
    }
};

void PrintTo(const QN& qn, std::ostream* os) {
    *os << "[n_x, n_y, n_z] = [ " << qn.n_x << ", " << qn.n_y << ", " << qn.n_z << "]";
}

std::vector<QN> generate_tests(std::size_t N) {
    std::vector<QN> tests;

    tests.resize(N*N*N);
    std::size_t index = 0;
    for (std::size_t i = 0; i < N; i++) {
        for (std::size_t j = 0; j < N; j++) {
            for (std::size_t k = 0; k < N; k++) {
                tests[index] = QN(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k));
                index++;
            }
        }
    }
    return tests;
}

class HarmonicOscillatorTest : public ::testing::TestWithParam<QN> {};

TEST_P(HarmonicOscillatorTest, QuantumNumbers) {
    HarmonicOscillator ho(1.0);

    QN qn = GetParam();
    std::size_t index = ho.get_index_from_quantum_numbers(qn.to_vec());
    ho.set_energy_level_index(index);
    EXPECT_STREQ(ho.get_state_string().c_str(), qn.to_str(index).c_str());
}

#ifdef FULL_TEST
INSTANTIATE_TEST_SUITE_P(SampleValues, HarmonicOscillatorTest, ::testing::ValuesIn(generate_tests(10)));
#else
INSTANTIATE_TEST_SUITE_P(SampleValues, HarmonicOscillatorTest, ::testing::ValuesIn(generate_tests(3)));
#endif

