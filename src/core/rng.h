#ifndef RNG_H
#define RNG_H

#include <random>

class RNG {
public:
    std::mt19937 m_engine;

    RNG();

    double uniform();
    double normal();

private:
    std::uniform_real_distribution<double> m_uniform;
    std::normal_distribution<double> m_normal;
};

#endif
