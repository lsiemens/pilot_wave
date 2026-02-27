#include "quantum/harmonic_oscillator.h"

#include <iostream>
#include <algorithm>
#include <iterator>

#include "core/math_util.h"

HarmonicOscillator::QuantumNumbers::QuantumNumbers() {
    m_n_x = 0;
    m_n_y = 0;
    m_n_z = 0;
}

HarmonicOscillator::QuantumNumbers::QuantumNumbers(std::size_t n_x, std::size_t n_y, std::size_t n_z) {
    m_n_x = n_x;
    m_n_y = n_y;
    m_n_z = n_z;
}

bool HarmonicOscillator::QuantumNumbers::isValid() const {
    return true;
}

bool HarmonicOscillator::QuantumNumbers::operator==(const QuantumNumbers& other) const {
    if (m_n_x != other.m_n_x) {
        return false;
    }

    if (m_n_y != other.m_n_y) {
        return false;
    }

    if (m_n_z != other.m_n_z) {
        return false;
    }

    return true;
}


HarmonicOscillator::HarmonicOscillator(double omega) : m_omega(omega), m_log_omega(std::log(omega)), m_sqrt_m_omega_hbar(std::sqrt(m_m_e*omega/m_hbar)) {
    // TODO solve for the widdth
    double N = 3; // s = 0, p = 1, d = 2, f = 3
    double r_rms = std::sqrt((2*N + 3)/(2*m_omega));
    m_width = 2*2*r_rms;
    m_origin = glm::dvec3(0., 0., 0.);
    
    m_norm = std::pow(m_m_e*m_omega/(PI_D*m_hbar), 0.75);

    find_energy_levels();
}

std::size_t HarmonicOscillator::get_index_from_quantum_numbers(std::vector<int> qn) {
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

std::string HarmonicOscillator::get_state_string() const {
    if (get_num_states() > 0) {
        std::string str_repr = "State magnatude [" + std::to_string(get_norm()) + "]";
        return str_repr;
    } else {
        std::size_t energy_level_index = get_energy_level_index();
        QuantumNumbers quantum_numbers = m_energy_levels_QN[energy_level_index];
        std::string str_repr = "Energy level [" + std::to_string(energy_level_index) + "] ";
        str_repr += "quantum numbers: ";
        str_repr += "(" + std::to_string(quantum_numbers.m_n_x) + ","
                        + std::to_string(quantum_numbers.m_n_y) + ","
                        + std::to_string(quantum_numbers.m_n_z) + ")";
        return str_repr;
    }
}

double HarmonicOscillator::psi_0_max() const {
    return m_norm;
}

double HarmonicOscillator::psi_n(glm::dvec3 position, std::size_t energy_level_index) const {
    QuantumNumbers quanum_numbers = m_energy_levels_QN[energy_level_index];

    return ( psi_n_1D(quanum_numbers.m_n_x, position.x)
            *psi_n_1D(quanum_numbers.m_n_y, position.y)
            *psi_n_1D(quanum_numbers.m_n_z, position.z));
}

glm::dvec3 HarmonicOscillator::grad_psi_n(glm::dvec3 position, std::size_t energy_level_index) const {
    QuantumNumbers quanum_numbers = m_energy_levels_QN[energy_level_index];

    double grad_x =( dpsi_n_1D_dx(quanum_numbers.m_n_x, position.x)
                    * psi_n_1D(quanum_numbers.m_n_y, position.y)
                    * psi_n_1D(quanum_numbers.m_n_z, position.z));

    double grad_y =(  psi_n_1D(quanum_numbers.m_n_x, position.x)
                    *dpsi_n_1D_dx(quanum_numbers.m_n_y, position.y)
                    * psi_n_1D(quanum_numbers.m_n_z, position.z));

    double grad_z =(  psi_n_1D(quanum_numbers.m_n_x, position.x)
                    * psi_n_1D(quanum_numbers.m_n_y, position.y)
                    *dpsi_n_1D_dx(quanum_numbers.m_n_z, position.z));

    return glm::dvec3(grad_x, grad_y, grad_z);
}

void HarmonicOscillator::find_energy_levels() {
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
}

double HarmonicOscillator::get_energy_eigenvalue(QuantumNumbers quantum_numbers) const {
    double N = static_cast<double>(quantum_numbers.m_n_x + quantum_numbers.m_n_y + quantum_numbers.m_n_z);

    return m_hbar*m_omega*(N + 1.5);
}

double HarmonicOscillator::hermite_n(std::size_t n, double x) const {
    if (n == 0) {
        return 1.;
    } else if (n == 1) {
        return 2*x;
    }

    double Hn, Hn_m1, Hn_m2;

    Hn_m1 = 1.;
    Hn = 2*x;
    for (std::size_t i = 1; i < n; i++) {
        Hn_m2 = Hn_m1;
        Hn_m1 = Hn;

        Hn = 2*(x*Hn_m1 - static_cast<double>(i)*Hn_m2);
    }

    return Hn;
}

double HarmonicOscillator::psi_n_1D(std::size_t n, double x) const {
    static const double log2 = std::log(2);
    static const double logPI = std::log(PI_D);

    double ln_prefactor = 0.25*(m_log_omega - logPI) - 0.5*(static_cast<double>(n)*log2 + ln_factorial(n));
    double exponent = ln_prefactor - (m_m_e*m_omega*x*x/(2.*m_hbar));
    return std::exp(exponent)*hermite_n(n, m_sqrt_m_omega_hbar*x);
}

double HarmonicOscillator::dpsi_n_1D_dx(std::size_t n, double x) const {
    static const double log2 = std::log(2);
    static const double logPI = std::log(PI_D);

    double ln_prefactor = 0.25*(m_log_omega - logPI) - 0.5*(static_cast<double>(n)*log2 + ln_factorial(n));
    double exponent = ln_prefactor - (m_m_e*m_omega*x*x/(2.*m_hbar));
    double polynomial = ((m_m_e*m_omega/m_hbar)*x*hermite_n(n, m_sqrt_m_omega_hbar*x)
                         -m_sqrt_m_omega_hbar*hermite_n(n + 1, m_sqrt_m_omega_hbar*x));
    return std::exp(exponent)*polynomial;
}

