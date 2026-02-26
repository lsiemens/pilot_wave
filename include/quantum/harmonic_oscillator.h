/// Find solutions to the harmonic oscillator as a product of the solutions to
/// the 1D. Linear combinations of degenerate solutions can be used to
/// reconstruct solutions in the basis of angular momentum states.

#ifndef HARMONIC_OSCILLATOR_H
#define HARMONIC_OSCILLATOR_H

#include <vector>
#include <algorithm>

#include <glm/glm.hpp>

#include "quantum/qstate.h"

/// Impliments the QState interface for the simple harmonic oscillator.
/// @ingroup quantum_mechanics
class HarmonicOscillator : public QState {
public:
    struct QuantumNumbers {
        std::size_t m_n_x;
        std::size_t m_n_y;
        std::size_t m_n_z;

        QuantumNumbers();
        QuantumNumbers(std::size_t n_x, std::size_t n_y, std::size_t n_z);
        bool isValid() const;
    };

    HarmonicOscillator(double omega);

    double psi_0_max() const override;
    double psi_n(glm::dvec3 position, std::size_t energy_level_index) const override;
    glm::dvec3 grad_psi_n(glm::dvec3 position, std::size_t energy_level_index) const override;
    std::size_t get_index_from_quantum_numbers(std::vector<int> qn) override;
    std::string get_state_string() const override;
    void find_energy_levels() override;

    double energy_eigenvalue(QuantumNumbers quantum_numbers) const;
    double factorial(std::size_t n) const;
    double hermite_n(std::size_t n, double x) const;
    double psi_n_1D(std::size_t n, double x) const;
    double dpsi_n_1D_dx(std::size_t n, double x) const;

private:
    std::vector<QuantumNumbers> m_energy_levels_QN;    
    std::vector<double> m_cumsum_ln_n;
    std::size_t m_N_max; // max degeneracy level
    double m_norm;
    const double m_omega;
    const double m_log_omega;
    const double m_sqrt_m_omega_hbar;
};

#endif
