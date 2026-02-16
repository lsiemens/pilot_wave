/// @defgroup quantum_mechanics Quantum Mechanics
/// @brief Overview of time dependent quantum mechanics.
///
/// INTRODUCTION ...
///
/// @section overview Overview of notation and Equations
/// Here I will give a quick overview of quantum mechanics, presenting the
/// equations and the notation that will be used. Latter we will move onto the
/// probability current and the unit system.
///
/// @subsection schrodingers_equation Schrodinger's Equation
/// The general form of Schrodinger's Equation if given by
/// \f$i\hbar \frac{\partial}{\partial t} \psi(\vec{x}, t) = \hat{H} \psi(\vec{x}, t)\f$,
/// where \f$\psi(\vec{x}, t)\f$ is the wavefunction in the position basis and
/// \f$\hat{H}\f$ is the Hamiltonian operator. In the position basis the
/// Hamiltonian typically have the form \f$\hat{H} =
/// \frac{-\hbar^2}{2 m}\nabla^2 + V(\vec{x})\f$, where \f$V(\vec{x})\f$ is the
/// potential of the system. The energy eigenvalues and eigenstates associated
/// with the Hamiltonian are defined by the equation
/// \f$\hat{H}\psi_n(\vec{x}) = E_n\psi_n(\vec{x})\f$, where \f$n\f$ is the
/// index of a state (ordered by the energy level), \f$E_n\f$ is the energy of
/// the \f$n\f$th energy level and \f$\psi_n(\vec{x})\f$ is the wave function of
/// the \f$n\f$th energy level. The time dependent solutions so Schrodinger's
/// equation can be expressed as
/// \f$\psi(\vec{x}, t) = \sum_{n} c_n e^{-i E_n t/\hbar} \psi_n(\vec{x})\f$,
/// where \f$c_n\f$ are complex numbers expressing the amplitudes of the
/// different energy levels at \f$t = 0\f$. Provided that the states
/// \f$\psi_n(\vec{x})\f$ are normalized, then the normalization of the
/// time dependent wavefunction is given by the condition
/// \f$\sum_n c_n^*c_n = 1\f$.
///
/// @subsection probability_current Probability Current
///
/// @section flow Flow
///
/// @section hartree_units Hartree Units
/// Calculations and equations are performed in Hartree atomic units, where
/// \f$\hbar = 1\f$,  \f$m_e = 1\f$, \f$e = 1\f$ and \f$4\pi\epsilon_0 = 1\f$.


#ifndef QSTATE_H
#define QSTATE_H

#include <vector>
#include <complex>

#include <glm/glm.hpp>

/// Public interface for quantum systems
/// @ingroup quantum_mechanics
class QState {
public:
    // --- Physical constants ---

    /// Reduced Planks's constant.
    /// The value of Plank's reduced constant is \f$\hbar = 1\f$ in units of
    /// \f$[\hbar]\f$.
    static constexpr double m_hbar = 1.;

    /// Mass of the electron.
    /// The mass of the electron is \f$m_e = 1\f$ in units of \f$[m_e]\f$.
    static constexpr double m_m_e = 1.;

    // --- Public methods

    virtual ~QState()=default;

    /// Set the initial state coefficients.
    /// Switch the mode to represent time dependent wave functions, the
    /// wavefunction is \f$\psi(\vec{x}, t) = \sum_{n}c_n e^{-i E_n t / \hbar}\psi_n(\vec{x})\f$.
    /// @param coefficients The initial coefficients \f$c_n\f$.
    void set_coefficients(const std::vector<std::complex<double>>& coefficients);

    /// Set the current energy level.
    /// Switch the mode to represent eigenstates, the wavefunction is
    /// \f$\psi_n(x)\f$.
    /// @param energy_level The index \f$n\f$ of the desired eigenstate.
    void set_energy_level(std::size_t energy_level);

    /// Update the quantum state.
    /// Increment the internal time and update the time dependent coefficient
    /// vectors coeff_t.
    /// @param dt The time step of the current update in units of
    /// \f$[\hbar / E_h]\f$.
    void update(double dt);

    /// Validate the state vectors.
    bool validate() const;

    // --- Primary physics output ---

    /// Characteristic probability density.
    /// Get the maximum probability density of the ground state. This density
    /// is given in units \f$[1 / {a_0}^3]\f$.
    double probability_density_0_max() const;

    /// The probability density scalar field.
    /// Evaluates the probability density
    /// \f$\rho(\vec{x}) = \psi^*(\vec{x}, t)\psi(\vec{x}, t)\f$.
    /// @param position The location \f$\vec{x}\f$ at which to evaluate the
    /// probability density, position is a 3 vector with units of \f$[a_0]\f$.
    /// @returns The probability density is given in units of
    /// \f$[1 / {a_0}^3]\f$.
    double probability_density(glm::dvec3 position) const;

    /// The probability current vector field.
    /// Evaluates the probability current
    /// \f$\vec{j}(\vec{x}, t) = \frac{\hbar}{m_e} \mathcal{J}\left( \psi^*(\vec{x}, t)\vec{\nabla}\psi(\vec{x}, t) \right)\f$.
    /// @param position The location \f$\vec{x}\f$ at which to evaluate the
    /// probability current, position is a 3 vector with units of \f$[a_0]\f$.
    /// @returns The probability current is given in units of
    /// \f$[a_0 E_h / \hbar]\f$.
    glm::dvec3 probability_current(glm::dvec3 position) const;

    // --- Virtual functions ---

    /// Characteristic probability amplitude.
    /// Get the maximum probability amplitude of the ground state. This
    /// amplitude is given in units \f$[{a_0}^{-3/2}]\f$.
    virtual double psi_0_max() const = 0;

    /// The Nth energy eigenstate.
    /// Evaluates the wavefunction \f$\psi_n(\vec{x})\f$.
    /// @param position The location \f$\vec{x}\f$ at which to evaluate the
    /// wavefunction, position is a 3 vector with units of \f$[a_0]\f$.
    /// @param energy_level The index \f$n\f$ of the energy eigenstate to
    /// evaluate.
    /// @returns The probability density of the state is given in units of
    /// \f$[{a_0}^{-3/2}]\f$.
    virtual double psi_n(glm::dvec3 position, std::size_t energy_level) const = 0;

    /// The gradient of the Nth energy eigenstate.
    /// Evaluates the gradient of the function \f$\vec{\nabla}\psi_n(\vec{x})\f$.
    /// @param position The location \f$\vec{x}\f$ at which to evaluate the
    /// gradient, position is a 3 vector with units of \f$[a_0]\f$.
    /// @param energy_level The index \f$n\f$ of the energy eigenstate to
    /// evaluate.
    /// @returns The gradient of the probability density of the state is given
    /// in units of \f$[{a_0}^{-5/2}]\f$.
    virtual glm::dvec3 grad_psi_n(glm::dvec3 position, std::size_t energy_level) const = 0;

    /// Compute the energy levels of the system.
    /// Find energy levels and the set of quantum numbers ordered by their
    /// energy. This is used to access energy eigenstates by their index.
    virtual void find_energy_levels() = 0;

    // --- Getters ---

    /// Get the number of coefficients.
    /// Get the number of states defined in the private vector m_coeff. If this
    /// is zero, then the system is in an energy eigenstates defined by the
    /// private index m_energy_level.
    std::size_t get_num_states() const {
        return m_num_states;
    }

    /// Get the energy level.
    /// Get the energy level of the system defined by the private index
    /// m_energy_level. The system only uses this index if the number of states
    /// in the private vector m_coeff is zero.
    std::size_t get_energy_level() const {
        return m_energy_level;
    }

    /// Get the width of the system.
    /// Typically this is an estimate based on the characteristic width of
    /// states such that there is a finite domain that contains most of the
    /// structure. The width is given in units \f$[a_0]\f$.
    double get_width() const {
        return m_width;
    }

    /// Get the origin of the system.
    /// Typically this is an estimate based on the mean position of states such
    /// that the finite domain contains most of the structure. The width is a 3
    /// vector given in units \f$[a_0]\f$.
    glm::dvec3 get_origin() const {
        return m_origin;
    }

protected:
    QState() = default;

    /// domain width
    double m_width = 0.;

    /// domain center
    glm::dvec3 m_origin = glm::dvec3(0., 0., 0.);

private:
    /// Elapsed time.
    /// The time used in the time dependent component of the wavefunction. The
    /// time is given in units of \f$[\hbar / E_h]\f$.
    double m_time = 0.;

    /// The systems energy levels.
    /// An ordered list of energy levels for use in the time dependent component
    /// of the wavefunction. The energy is given in units of \f$[E_h]\f$.
    std::vector<double> m_energy_eigenvalues = {};

    /// The initial state coefficients.
    /// A list of the initial complex coefficients ordered by the energy level.
    std::vector<std::complex<double>> m_coeff = {};

    /// The time dependant state coefficients.
    /// A list of the time dependant complex coefficients orderd by the energy
    /// level.
    std::vector<std::complex<double>> m_coeff_t = {};

    /// The energy level of the current energy eigenstate.
    /// The index of the energy level for the system. Used if m_num_states is
    /// zero.
    std::size_t m_energy_level = 0;

    /// The number of defined state coefficients.
    /// This is given by m_coeff.size() and determines if the system is
    /// representing a time depended state or one of the energy eigenstates.
    std::size_t m_num_states = 0;
};

#endif
