#include "quantum/hydrogen_atom.h"

#include <iostream>
#include <algorithm>
#include <iterator>

#include "core/math_util.h"

HydrogenAtom::QuantumNumbers::QuantumNumbers() {
    m_n = 1;
    m_l = 0;
    m_m = 0;
}

HydrogenAtom::QuantumNumbers::QuantumNumbers(std::size_t n, std::size_t l, int m) {
    m_n = n;
    m_l = l;
    m_m = m;
}

bool HydrogenAtom::QuantumNumbers::isValid() const {
    if (m_n == 0) {
        return false;
    }

    if (m_l >= m_n) {
        return false;
    }

    if (static_cast<std::size_t>(std::abs(m_m)) > m_l) {
        return false;
    }

    return true;
}

bool HydrogenAtom::QuantumNumbers::operator==(const QuantumNumbers& other) const {
    if (m_n != other.m_n) {
        return false;
    }

    if (m_l != other.m_l) {
        return false;
    }

    if (m_m != other.m_m) {
        return false;
    }

    return true;
}


HydrogenAtom::HydrogenAtom(std::size_t Z, double mu) {
    find_energy_levels();
}

std::size_t HydrogenAtom::get_index_from_quantum_numbers(std::vector<int> qn) {
    if (qn.size() < 3) {
        throw std::range_error("Three quantum numbers are required.");
        // raise exception
    }

    QuantumNumbers quantum_numbers(qn[0], qn[1], qn[2]);

    if (not quantum_numbers.isValid()) {
        throw std::runtime_error("Invalid set of quantum numbers.");
    }

    std::size_t num_states = get_num_states();
    double energy_level = get_energy_eigenvalue(quantum_numbers);
    if (energy_level < m_energy_eigenvalues[num_states - 1]) {
        for (std::size_t i = 0; i < num_states; i++) {
            if (quantum_numbers == m_energy_levels_QN[i]) {
                return i;
            }
        }

        throw std::runtime_error("There are missing energy levels.");
    } else {
        set_coefficient(2*num_states, {0., 0.});
        return get_index_from_quantum_numbers(qn);
    }
}

std::string HydrogenAtom::get_state_string() const {
    if (get_num_states() > 0) {
        std::string str_repr = "State magnatude [" + std::to_string(get_norm()) + "]";
        return str_repr;
    } else {
        std::size_t energy_level_index = get_energy_level_index();
        QuantumNumbers quantum_numbers = m_energy_levels_QN[energy_level_index];
        std::string str_repr = "Energy level [" + std::to_string(energy_level_index) + "] ";
        str_repr += "quantum numbers: ";
        str_repr += "(" + std::to_string(quantum_numbers.m_n) + ","
                        + std::to_string(quantum_numbers.m_l) + ","
                        + std::to_string(quantum_numbers.m_m) + ")";
        return str_repr;
    }
}

double HydrogenAtom::psi_0_max() const {
    return m_norm;
}

double HydrogenAtom::psi_n(glm::dvec3 position, std::size_t energy_level_index) const {
    QuantumNumbers quanum_numbers = m_energy_levels_QN[energy_level_index];

    return 0.;
}

glm::dvec3 HydrogenAtom::grad_psi_n(glm::dvec3 position, std::size_t energy_level_index) const {
    QuantumNumbers quanum_numbers = m_energy_levels_QN[energy_level_index];

    return glm::dvec3(0., 0., 0.);
}

void HydrogenAtom::find_energy_levels() {
    /*
    std::size_t n_states = 0;
    if (get_num_states() > 0) {
        n_states = get_num_states();
    } else {
        n_states = get_energy_level_index() + 1;
    }

    //std::cout << "Find energy levels: # states " << n_states << std::endl;

    std::size_t N_max = static_cast<std::size_t>(std::ceil(std::pow(6.*static_cast<double>(n_states), 1/3.) - 2.));

    while ((N_max + 1)*(N_max + 2)*(N_max + 3)/6 < n_states) {
        N_max++;
    }

    std::vector<QuantumNumbers> initial_set((N_max + 1)*(N_max + 2)*(N_max + 3)/6);

    std::size_t index = 0;
    for (std::size_t n_z = 0; n_z <= N_max; n_z++) {
        for (std::size_t n_y = 0; n_y <= N_max - n_z; n_y++) {
            for (std::size_t n_x = 0; n_x <= N_max - n_z - n_y; n_x++) {
                initial_set[index] = QuantumNumbers(n_x, n_y, n_z);
                index++;
            }
        }
    }

    // TODO sort by n^2 then by E using a stable sorting algorithm so degenerate
    // blocks are grouped by spin
    auto comparison = [this](const QuantumNumbers& a, const QuantumNumbers& b) {
        double E_a = get_energy_eigenvalue(a);
        double E_b = get_energy_eigenvalue(b);

        if (E_a == E_b) {
            if (a.m_n_x != b.m_n_x) {
                return a.m_n_x < b.m_n_x;
            }

            if (a.m_n_y != b.m_n_y) {
                return a.m_n_y < b.m_n_y;
            }

            if (a.m_n_z != b.m_n_z) {
                return a.m_n_z < b.m_n_z;
            }

            throw std::runtime_error("No ordering, quanum numbers are equal");
        }

        return E_a < E_b;
    };
    std::sort(initial_set.begin(), initial_set.end(), comparison);

    m_energy_levels_QN.resize(n_states);
    std::copy_n(initial_set.begin(), n_states, m_energy_levels_QN.begin());

    m_energy_eigenvalues.resize(n_states);
    for (std::size_t i = 0; i < n_states; i++) {
        m_energy_eigenvalues[i] = get_energy_eigenvalue(m_energy_levels_QN[i]);
    }

    m_N_max = N_max;
    */
}

double HydrogenAtom::get_energy_eigenvalue(QuantumNumbers quantum_numbers) const {
    double N = static_cast<double>(quantum_numbers.m_n);

    return 0.;
}
