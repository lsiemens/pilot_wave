#include "quantum/square_well.h"

#include "testing/testing.h"
#include <gtest/gtest.h>

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
    for (std::size_t i = 1; i < N + 1; i++) {
        for (std::size_t j = 1; j < N + 1; j++) {
            for (std::size_t k = 1; k < N + 1; k++) {
                tests[index] = QN(static_cast<int>(i), static_cast<int>(j), static_cast<int>(k));
                index++;
            }
        }
    }
    return tests;
}

class SquareWellTest : public ::testing::TestWithParam<QN> {};

TEST_P(SquareWellTest, QuantumNumbers) {
    SquareWell sw(1.0);

    QN qn = GetParam();
    std::size_t index = sw.level_from_quantum_numbers(qn.to_vec());
    sw.set_energy_level(index);
    EXPECT_STREQ(sw.get_state_string().c_str(), qn.to_str(index).c_str());
}

#ifdef FULL_TEST
INSTANTIATE_TEST_SUITE_P(SampleValues, SquareWellTest, ::testing::ValuesIn(generate_tests(10)));
#else
INSTANTIATE_TEST_SUITE_P(SampleValues, SquareWellTest, ::testing::ValuesIn(generate_tests(3)));
#endif
