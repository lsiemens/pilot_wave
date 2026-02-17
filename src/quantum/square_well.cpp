#include "quantum/square_well.h"

#include <iostream>

#include "core/math_util.h"

double dsin3_dx(double kx, double x, double ky, double y, double kz, double z) {
    return kx*std::cos(kx*x)*std::sin(ky*y)*std::sin(kz*z);
}

SquareWell::QuantumNumbers::QuantumNumbers() {
    m_n_x = 1;
    m_n_y = 1;
    m_n_z = 1;
}

SquareWell::QuantumNumbers::QuantumNumbers(std::size_t n_x, std::size_t n_y, std::size_t n_z) {
    m_n_x = n_x;
    m_n_y = n_y;
    m_n_z = n_z;
}

bool SquareWell::QuantumNumbers::isValid() const {
    return ((m_n_x > 0) and (m_n_y > 0) and (m_n_z > 0));
}

SquareWell::SquareWell(double width) {
    m_width = width;
    m_origin = glm::dvec3(width/2., width/2., width/2.);
    
    m_norm = std::sqrt(2./m_width)*2./m_width;

    find_energy_levels();
}

double SquareWell::psi_0_max() const {
    return m_norm;
}

double SquareWell::psi_n(glm::dvec3 position, std::size_t energy_level) const {
    assert(energy_level < energy_levels_QN.size());
    QuantumNumbers quantum_numbers = energy_levels_QN[energy_level];

    if ((position.x < 0) or (position.x > m_width)) {
        return 0.;
    }

    if ((position.y < 0) or (position.y > m_width)) {
        return 0.;
    }

    if ((position.z < 0) or (position.z > m_width)) {
        return 0.;
    }

    double k_x = PI_D*static_cast<double>(quantum_numbers.m_n_x)/m_width;
    double k_y = PI_D*static_cast<double>(quantum_numbers.m_n_y)/m_width;
    double k_z = PI_D*static_cast<double>(quantum_numbers.m_n_z)/m_width;

    return m_norm*std::sin(k_x*position.x)
                 *std::sin(k_y*position.y)
                 *std::sin(k_z*position.z);
}

glm::dvec3 SquareWell::grad_psi_n(glm::dvec3 position, std::size_t energy_level) const {
    assert(energy_level < energy_levels_QN.size());
    QuantumNumbers quantum_numbers = energy_levels_QN[energy_level];

    if ((position.x < 0) or (position.x > m_width)) {
        return glm::dvec3(0., 0., 0.);
    }

    if ((position.y < 0) or (position.y > m_width)) {
        return glm::dvec3(0., 0., 0.);
    }

    if ((position.z < 0) or (position.z > m_width)) {
        return glm::dvec3(0., 0., 0.);
    }

    double k_x = PI_D*static_cast<double>(quantum_numbers.m_n_x)/m_width;
    double k_y = PI_D*static_cast<double>(quantum_numbers.m_n_y)/m_width;
    double k_z = PI_D*static_cast<double>(quantum_numbers.m_n_z)/m_width;

    glm::dvec3 grad_psi;
    grad_psi.x = m_norm*dsin3_dx(k_x, position.x, k_y, position.y, k_z, position.z);
    grad_psi.y = m_norm*dsin3_dx(k_y, position.y, k_z, position.z, k_x, position.x);
    grad_psi.z = m_norm*dsin3_dx(k_z, position.z, k_x, position.x, k_y, position.y);
    return grad_psi;
}

void SquareWell::find_energy_levels() {
    std::size_t n_states = 0;
    if (get_num_states() > 0) {
        n_states = get_num_states();
    } else {
        n_states = get_energy_level() + 1;
    }

    //std::cout << "Find energy levels: # states " << n_states << std::endl;

    // TODO improve and clarify this huristic
    std::size_t n_max = static_cast<std::size_t>(std::sqrt(3.) + std::ceil(std::pow(6.*static_cast<double>(n_states)/3.14, 1/3.)));
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

    energy_levels_QN.resize(n_states);
    std::copy_n(initial_set.begin(), n_states, energy_levels_QN.begin());
}

double SquareWell::energy_eigenvalue(QuantumNumbers quantum_numbers) const {
    if (not quantum_numbers.isValid()) {
        throw std::out_of_range("Invalid quantum numbers");
    }

    double k_x = PI_D*static_cast<double>(quantum_numbers.m_n_x)/m_width;
    double k_y = PI_D*static_cast<double>(quantum_numbers.m_n_y)/m_width;
    double k_z = PI_D*static_cast<double>(quantum_numbers.m_n_z)/m_width;

    double k_sqr = k_x*k_x + k_y*k_y + k_z*k_z;
    return m_hbar*m_hbar*k_sqr/(2.*m_m_e);
}
