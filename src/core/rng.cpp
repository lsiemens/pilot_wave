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

/// Sample from Poisson distribution using inverse sampling
std::size_t RNG::poisson(double lambda) {
    std::size_t n = 0;
    double p_n = std::exp(-lambda);
    double cumsum_p_n = p_n;
    double sample = uniform();
    
    while (sample > cumsum_p_n) {
        n++;
        p_n *= lambda/static_cast<double>(n);
        cumsum_p_n += p_n;
    }
    return n;
}
