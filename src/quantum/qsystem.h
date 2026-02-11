#ifndef QSYSTEM_H
#define QSYSTEM_H

#include <limits>
#include <array>
#include <concepts>
#include <complex>
#include <stdexcept>

#include <glm/glm.hpp>

// Only here temeraily for testing
using namespace std::complex_literals;

template<std::size_t N_STATES>
using Coefficients = std::array<std::complex<double>, N_STATES>;

template <typename QN>
concept IsQuantumNumbers = 
    requires(const QN& quantum_numbers) {
        { quantum_numbers.isValid() } -> std::same_as<bool>;
    };

template <typename P>
using QN_type = typename P::QuantumNumbers;

template <typename P, std::size_t N_STATES>
concept IsPotential = 
    requires(P potential) {
        typename P::QuantumNumbers;
    }
    and IsQuantumNumbers<QN_type<P>>
    and requires(const P& potential, glm::vec3 position, QN_type<P> quantum_numbers) {
        { potential.V(position) } -> std::same_as<double>;
        { potential.energy_eigenvalue(quantum_numbers) } -> std::same_as<double>;
        { potential.energy_eigenstate(quantum_numbers, position) } -> std::same_as<double>;
        { potential.template get_levels<N_STATES>() } -> std::same_as<std::array<QN_type<P>, N_STATES>>;
    };

template <typename Potential, std::size_t N_STATES>
requires IsPotential<Potential, N_STATES>
class QSystem {
public:
    using QuantumNumbers = QN_type<Potential>;

    const Potential m_potential;
    const std::array<QuantumNumbers, N_STATES> m_levels;

    QSystem(const Potential& potential) : m_potential(potential), m_levels(potential.template get_levels<N_STATES>()) {}

    double state(std::size_t n, glm::vec3 position) const {
        if (n >= N_STATES) {
            throw std::out_of_range("State (n = " + std::to_string(n) + ") out of range [0, N_STATES - 1] with N_STATES = " + std::to_string(N_STATES) + ".");
        }

        return m_potential.energy_eigenstate(m_levels[n], position);
    }

    double energy(std::size_t n) const {
        if (n >= N_STATES) {
            throw std::out_of_range("State (n = " + std::to_string(n) + ") out of range [0, N_STATES - 1] with N_STATES = " + std::to_string(N_STATES) + ".");
        }

        return m_potential.energy_eigenvalue(m_levels[n]);
    }
};

#endif
