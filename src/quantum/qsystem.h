#ifndef QSYSTEM_H
#define QSYSTEM_H

#include <limits>
#include <vector>
#include <concepts>
#include <complex>

#include <glm/glm.hpp>

using Coefficients = std::vector<std::complex<double>>;

// Only here temeraily for testing
using namespace std::complex_literals;

template <typename QN>
concept IsQuantumNumbers = 
    requires(const QN& quantum_numbers) {
        { quantum_numbers.isValid() } -> std::same_as<bool>;
    };

template <typename P>
using QN_type = typename P::QuantumNumbers;

template <typename P>
concept IsPotential = 
    requires(P potential) {
        typename P::QuantumNumbers;
    }
    and IsQuantumNumbers<QN_type<P>>
    and requires(const P& potential, glm::vec3 position, QN_type<P> quantum_numbers) {
        { potential.V(position) } -> std::same_as<double>;
        { potential.energy_eigenvalue(quantum_numbers) } -> std::same_as<double>;
        { potential.energy_eigenstate(quantum_numbers, position) } -> std::same_as<double>;
    };


template <IsPotential Potential>
class QSystem {
public:
    using QuantumNumbers = QN_type<Potential>;

    const Potential m_potential;
    Coefficients m_coefficients;

    QSystem(const Potential& potential, const Coefficients& coefficients) : m_potential(potential), m_coefficients(coefficients) {}

    double wave_function(glm::vec3 position) const {
        double value = 0.;
        for (std::size_t i = 0; i < m_coefficients.size(); i++) {
            QuantumNumbers qn = QuantumNumbers(i + 1);
            value += m_coefficients[i]*m_potential.energy_eigenstate(qn, position);
        }
        return value;
    }
};


class MyPotential {
public:
    struct QuantumNumbers {
        int n;
        constexpr bool isValid() const {
            return (n > 0);
        }
    };

    constexpr double V(glm::vec3 position) const {
        if ((0. <= position.x) and (position.x < 1.)) {
            return 0.;
        } else {
            return std::numeric_limits<double>::infinity();
        }
    }

    constexpr double energy_eigenvalue(QuantumNumbers quantum_numbers) const {
        assert(quantum_numbers.isValid());

        return 1.0*static_cast<double>(quantum_numbers.n*quantum_numbers.n);
    }

    constexpr double energy_eigenstate(QuantumNumbers quantum_numbers, glm::vec3 position) const {
        assert(quantum_numbers.isValid());

        return std::sin(static_cast<double>(quantum_numbers.n)*position.x);
    }
};

// only here temperaily for test
QSystem<MyPotential> A(MyPotential(), {1. + 0.i, 0. + 2.i, 3., 4.i});

#endif
