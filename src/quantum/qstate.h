/// Using hartree units, in this system hbar = e = m_e = permitivity = 1

#ifndef QSTATE_H
#define QSTATE_H

#include <vector>
#include <complex>

#include <glm/glm.hpp>

class QState {
public:
    static constexpr double m_hbar = 1.;
    static constexpr double m_m_e = 1.;
    virtual ~QState()=default;

    void set_coefficients(const std::vector<std::complex<double>>& coefficients);
    void set_energy_level(std::size_t energy_level);
    void update(double dt);
    bool validate() const;

    double probability_density_0_max() const;
    double probability_density(glm::dvec3 position) const;
    glm::dvec3 probability_current(glm::dvec3 position) const;

    virtual double psi_0_max() const = 0;
    virtual double psi_n(glm::dvec3 position, std::size_t energy_level) const = 0;
    virtual glm::dvec3 grad_psi_n(glm::dvec3 position, std::size_t energy_level) const = 0;
    virtual void find_energy_levels() = 0;


    std::size_t get_num_states() const {
        return m_num_states;
    }

    std::size_t get_energy_level() const {
        return m_energy_level;
    }

    double get_width() const {
        return m_width;
    }

    glm::dvec3 get_origin() const {
        return m_origin;
    }

protected:
    QState() = default;

    double m_width;
    glm::dvec3 m_origin;

private:
    double m_time = 0.;
    std::vector<double> m_energy_eigenvalues;
    std::vector<std::complex<double>> m_coeff;
    std::vector<std::complex<double>> m_coeff_t;
    std::size_t m_energy_level;
    std::size_t m_num_states;
};

#endif
