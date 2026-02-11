#ifndef SQUARE_WELL_H
#define SQUARE_WELL_H

#include <vector>
#include <algorithm>
#include <complex>

#include <iostream>

#include "quantum/qsystem.h"

class MyPotential {
public:
    struct QuantumNumbers {
        int n_x = 1;
        int n_y = 1;
        int n_z = 1;

        constexpr bool isValid() const {
            return ((n_x > 0) and (n_y > 0) and (n_y > 0));
        }
    };

    constexpr double V(glm::vec3 position) const {
        if (    (0. <= position.x) and (position.x < 1.)
            and (0. <= position.y) and (position.y < 1.)
            and (0. <= position.z) and (position.z < 1.)) {
                return 0.;
        } else {
            return std::numeric_limits<double>::infinity();
        }
    }

    constexpr double energy_eigenvalue(QuantumNumbers quantum_numbers) const {
        if (not quantum_numbers.isValid()) {
            throw std::out_of_range("Invalid quantum numbers");
        }

        double n_x = static_cast<double>(quantum_numbers.n_x);
        double n_y = static_cast<double>(quantum_numbers.n_y);
        double n_z = static_cast<double>(quantum_numbers.n_z);
        
        return n_x*n_x + n_y*n_y + n_z*n_z;
    }

    constexpr double energy_eigenstate(QuantumNumbers quantum_numbers, glm::vec3 position) const {
        if (not quantum_numbers.isValid()) {
            throw std::out_of_range("Invalid quantum numbers");
        }

        double n_x = static_cast<double>(quantum_numbers.n_x);
        double n_y = static_cast<double>(quantum_numbers.n_y);
        double n_z = static_cast<double>(quantum_numbers.n_z);
        
        return std::sin(3.14*(n_x*position.x + n_y*position.y + n_z*position.z));
    }

    template <std::size_t N_STATES>
    std::array<QuantumNumbers, N_STATES> get_levels() const {
        std::array<QuantumNumbers, N_STATES> levels;

        std::size_t n_max = static_cast<std::size_t>(std::sqrt(3.) + std::ceil(std::pow(6.*N_STATES/3.14, 1/3.)));
        std::vector<QuantumNumbers> initial_set(n_max*n_max*n_max);

        for (std::size_t i = 0; i < n_max; i++) {
            for (std::size_t j = 0; j < n_max; j++) {
                for (std::size_t k = 0; k < n_max; k++) {
                    initial_set[k + n_max*j + n_max*n_max*i] = QuantumNumbers(i + 1, j + 1, k + 1);
                }
            }
        }

        auto comparison = [this](const QuantumNumbers& a, const QuantumNumbers& b) {
            return energy_eigenvalue(a) < energy_eigenvalue(b);
        };
        std::sort(initial_set.begin(), initial_set.end(), comparison);
        std::copy_n(initial_set.begin(), N_STATES, levels.begin());
        return levels;
    }
};

QSystem<MyPotential, 100> A = QSystem<MyPotential, 100>(MyPotential());

#endif
