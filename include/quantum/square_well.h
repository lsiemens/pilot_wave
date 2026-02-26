#ifndef SQUARE_WELL_H
#define SQUARE_WELL_H

#include <vector>
#include <algorithm>

#include <glm/glm.hpp>

#include "quantum/qstate.h"

/// Impliments the QState interface for the infinite square well potential.
/// @ingroup quantum_mechanics
class SquareWell : public QState {
public:
     struct QuantumNumbers {
        std::size_t m_n_x;
        std::size_t m_n_y;
        std::size_t m_n_z;

        QuantumNumbers();
        QuantumNumbers(std::size_t n_x, std::size_t n_y, std::size_t n_z);
        bool isValid() const;
        bool operator==(const QuantumNumbers& other) const;
    };

    SquareWell(double width);

    double psi_0_max() const override;
    double psi_n(glm::dvec3 position, std::size_t energy_level_index) const override;
    glm::dvec3 grad_psi_n(glm::dvec3 position, std::size_t energy_level_index) const override;
    void find_energy_levels() override;
    std::size_t get_index_from_quantum_numbers(std::vector<int> qn) override;
    std::string get_state_string() const override;
    double get_energy_eigenvalue(QuantumNumbers quantum_numbers) const;

private:
    std::vector<QuantumNumbers> m_energy_levels_QN;    
    double m_norm;
};

#endif
