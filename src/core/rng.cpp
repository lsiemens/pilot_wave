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
    static const double gaussian_cutoff = 500.;
    static const std::size_t max_itteration = static_cast<std::size_t>(gaussian_cutoff + 5*std::sqrt(gaussian_cutoff));
    std::size_t n = 0;

    if (lambda > gaussian_cutoff) {
        return static_cast<std::size_t>(lambda + m_normal(m_engine)*std::sqrt(lambda));
    } else {
        double p_n = std::exp(-lambda);
        double cumsum_p_n = p_n;
        double sample = uniform();

        while (sample > cumsum_p_n) {
            if (n >= max_itteration) {
                break;
            }

            n++;
            p_n *= lambda/static_cast<double>(n);
            cumsum_p_n += p_n;
        }
        return n;
    }
}
