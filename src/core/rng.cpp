#include "core/rng.h"

RNG::RNG() {
    std::random_device rd;
    std::seed_seq sequence{rd(), rd(), rd(), rd()};

    m_engine = std::mt19937(sequence);

    m_uniform = std::uniform_real_distribution<double>(0.0, 1.0);
    m_normal = std::normal_distribution<double>(0.0, 1.0);
}

double RNG::uniform() {
    return m_uniform(m_engine);
}

double RNG::normal() {
    return m_normal(m_engine);
}
