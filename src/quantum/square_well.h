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
        int n_x = 1;
        int n_y = 1;
        int n_z = 1;

        bool isValid() const {
            return ((n_x > 0) and (n_y > 0) and (n_z > 0));
        }
    };

    SquareWell(double width);

    double psi_0_max() const override;
    double psi_n(glm::dvec3 position, std::size_t energy_level) const override;
    glm::dvec3 grad_psi_n(glm::dvec3 position, std::size_t energy_level) const override;
    void find_energy_levels() override;

private:
    std::vector<QuantumNumbers> energy_levels_QN;    
    double m_norm;

    double energy_eigenvalue(QuantumNumbers quantum_numbers) const;
};

#endif
