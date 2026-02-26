#include "quantum/square_well.h"
#include "quantum/harmonic_oscillator.h"

#include "testing/testing.h"
#include "testing/integrate.h"
#include <gtest/gtest.h>

template <typename T>
class QStateTest : public ::testing::Test {
};

using QStateTypes = ::testing::Types<SquareWell, HarmonicOscillator>;
TYPED_TEST_SUITE(QStateTest, QStateTypes);

TYPED_TEST(QStateTest, Constants) {
    TypeParam state(1.0);
    state.set_energy_level(1);

    EXPECT_PRED_FORMAT2(IsClose, state.m_hbar, 1.0);
    EXPECT_PRED_FORMAT2(IsClose, state.m_m_e, 1.0);
}

TYPED_TEST(QStateTest, Orthonormal3D) {
    double param = 1.0;
#ifdef FULL_TEST
    std::size_t N_max = 10;
#else
    std::size_t N_max = 3;
#endif
    TypeParam state(param);
    state.set_energy_level(N_max);
    double x_max = 10.;

    std::size_t resolution = 100;
    for (std::size_t n = 0; n < N_max; n++) {
        for (std::size_t m = 0; m <= n; m++) {

            auto integrand_3D = [&state, n, m](glm::dvec3 x) {
                return state.psi_n(x, n)*state.psi_n(x, m);
            };

            auto value = integrate_uniform_3D(-x_max, x_max, resolution, integrand_3D);

            std::string name = "integral(-10, 10) |psi_3D(" + std::to_string(n) + ", x)psi_3D(" + std::to_string(m) + ", x)|^2";
            if (n == m) {
                EXPECT_PRED_FORMAT2(IsClose, value, 1.) << name;
            } else {
                EXPECT_PRED_FORMAT2(IsClose, value, 0.) << name;
            }
        }
    }
}
