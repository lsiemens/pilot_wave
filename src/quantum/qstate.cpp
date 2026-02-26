#include "quantum/qstate.h"

void QState::set_coefficients(const std::vector<std::complex<double>>& coefficients) {
    m_coeff = coefficients;
    m_coeff_t = coefficients;
    m_num_states = m_coeff.size();
    m_energy_level = 0;
    m_time = 0.;

    find_energy_levels();
    norm_update();
}

void QState::set_coefficient(std::size_t n, std::complex<double> coefficient) {
    if (n >= m_num_states) {
        m_coeff.resize(n + 1);
        for (std::size_t i = m_num_states; i < n + 1; i++) {
            m_coeff[i] = std::complex<double>(0., 0.);
        }
    }
    m_coeff[n] = coefficient;

    m_coeff_t = m_coeff;
    m_num_states = m_coeff.size();
    m_energy_level = 0;
    m_time = 0.;

    find_energy_levels();
    norm_update();
}

void QState::set_coefficient(std::vector<int> quantum_numbers, std::complex<double> coefficient) {
    std::size_t index = level_from_quantum_numbers(quantum_numbers);
    set_coefficient(index, coefficient);
}

void QState::normalize() {
    norm_update();
    double value = 1/std::sqrt(m_state_norm);
    for (std::size_t i = 0; i < m_num_states; i++) {
        m_coeff[i] = m_coeff[i]*value;
    }

    m_coeff_t = m_coeff;
    m_time = 0.;
    norm_update();
}

void QState::set_energy_level(std::size_t energy_level) {
    m_energy_level = energy_level;
    m_num_states = 0;

    find_energy_levels();
    norm_update();
}

void QState::update(double dt) {
    m_time += dt;
    
    for (std::size_t i = 0; i < m_num_states; i++) {
        double argument = m_energy_eigenvalues[i]*m_time/m_hbar;
        m_coeff_t[i] = m_coeff[i]*std::exp(std::complex<double>(0., -1.)*argument);
    }
}

double QState::probability_density_0_max() const {
    return psi_0_max()*psi_0_max();
}

double QState::probability_density(glm::dvec3 position) const {
    if (m_num_states > 0) {
        std::complex<double> psi = 0;
        for (std::size_t i = 0; i < m_num_states; i++) {
            psi += m_coeff_t[i]*psi_n(position, i);
        }
        return std::norm(psi);
    } else {
        double psi = psi_n(position, m_energy_level);
        return psi*psi;
    }
}

glm::dvec3 QState::probability_current(glm::dvec3 position) const {
    if (m_num_states > 0) {
        glm::dvec3 current = glm::dvec3(0.f, 0.f, 0.f);
        std::complex<double> psi_conj = 0;
        for (std::size_t i = 0; i < m_num_states; i++) {
            psi_conj += std::conj(m_coeff_t[i])*psi_n(position, i);
        }

        for (std::size_t i = 0; i < m_num_states; i++) {
            current += (m_hbar/m_m_e)*std::imag(psi_conj*m_coeff_t[i])*grad_psi_n(position, i);
        }

        return current;
    } else {
        return glm::dvec3(0.f, 0.f, 0.f);
    }
}

bool QState::validate() const {
    std::size_t num_states = m_coeff.size();

    if (m_energy_eigenvalues.size() != num_states) {
        return false;
    }

    if (m_coeff_t.size() != num_states) {
        return false;
    }

    return true;
}
