#include "core/math_util.h"

#include "testing/testing.h"
#include <gtest/gtest.h>

std::vector<std::size_t> generate_tests(std::size_t N) {
    std::vector<std::size_t> tests;

    tests.resize(N);
    for (std::size_t i = 0; i < N; i++) {
        tests[i] = i;
    }
    return tests;
}

class MathUtilTest : public ::testing::TestWithParam<std::size_t> {};

TEST_P(MathUtilTest, factorial) {
    std::size_t n = GetParam();

    double value = 1.;
    for (std::size_t i = 1; i <= n; i++) {
        value *= static_cast<double>(i);
    }
    std::string name = "factorial(" + std::to_string(n) + ")";
    EXPECT_PRED_FORMAT2(IsClose, factorial(n), value) << name;
}

TEST_P(MathUtilTest, ln_factorial) {
    std::size_t n = GetParam();

    double value = 1.;
    for (std::size_t i = 1; i <= n; i++) {
        value *= static_cast<double>(i);
    }
    std::string name = "ln(factorial(" + std::to_string(n) + "))";
    EXPECT_PRED_FORMAT2(IsClose, ln_factorial(n), std::log(value)) << name;
}

INSTANTIATE_TEST_SUITE_P(SampleValues, MathUtilTest, ::testing::ValuesIn(generate_tests(10)));
