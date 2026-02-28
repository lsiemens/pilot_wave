/// @defgroup hydrogen_atom The Hydrogen Atom
/// @brief Solution of the hydrogen atom in cartesian corrdinates
///
/// @section hydrogen_like_atom The Hydrogen Like Atom
/// The potential of a hydrogen like atom (an ionized atom such that it has one
/// electron) is described by the charge of the nucleaus \f$Z\f$ and is given by
/// \f$V(r) = -\frac{1}{4 \pi \epsilon_0}\frac{Ze^2}{r}\f$ where
/// \f$\epsilon_0\f$ is the permittivity of the vacuum. Classically this problem
/// has Keplerian orbits as solutions. In quantum mechanics shrodingers equation
/// with this potential is
/// \f$i\hbar\partial_t \psi(\vec(x)) = \frac{-\hbar^2}{2m}\nabla^2 \psi(\vec(x)) - \frac{1}{4\pi\epsilon_0}\frac{Ze^2}{r}\psi(\vec{x})\f$.
/// This problem has energy eigenfunctions that are seperated into a radial and
/// spherical part, such that energy eigenfunctions are \f$\psi_{n \ell m}(r, \theta, \phi) = R_{n \ell}(r)Y_{\ell m}(\theta, \phi)\f$,
/// where \f$Y_{\ell m}(\theta, \phi)\f$ are the spherical harmonics.
///
/// @section the_problem The Problem
/// The solution to the hydrogen states as given above is in spherical
/// coordinates and can be solved explicity interms of *generalized Laguerre
/// polynomials*, *associated Legendre polynomials* and some exponential terms.
/// These polynomials are functions of the radius and the cosine of \f$\theta\f$
/// respectivly. Computing these functions in cartesian coordinates involves
/// inverting trigometric functions which can be slow but the bigger issue is
/// the coordinate singularity at \f$r = 0\f$. This can be resolved by using
/// regular solid harmonics which are defined in cartesian coordinates and are
/// multivariate polinomials of \f$(x, y, z)\f$. The relation between the two
/// types of harmonics is
/// \f$\mathcal{Y}_{\ell m}(\vec{r}) = \sqrt{\frac{4\pi}{2\ell + 1}}r^\ell Y_{\ell m}(\theta, \phi)\f$
///
/// Solid harmonics are solutions of the differential equation
/// \f$\nabla^2 f(x, y, z) = 0\f$ while spherical harmonics are the solution to
/// the radial part of the equation \f$\nabla^2 f(r, \theta, \phi) = 0\f$. By
/// Seperation of variables the radial part is
/// \f$\frac{d}{dr}\left( r^2 \frac{d}{dr}R(r) \right) = \ell\left(\ell  + 1\right)R(r)\f$,
/// which has solutions \f$R(r) = Ar^\ell + Br^{-\ell - 1}\f$. So the solutions
/// to the solid harmonics are the regular solid harmonics given by
/// \f$r^\ell Y_{\ell m}(\theta, \phi)\f$ which are regular at \f$r = 0\f$ and
/// the irregular solid harmonics \f$r^{-\ell - 1}Y_{\ell m}(\theta, \phi)\f$ as
/// the name implies they are not regular at \f$r = 0\f$ and satisfy
/// \f$\nabla^2 f(r, \theta, \phi) = 0\f$ for \f$r \neq 0\f$. In cartesian
/// coordinates the regular solid harmonics are the functions that satisfy
/// \f$\nabla^2 f(x, y, z) = 0\f$, are regular at \f$r = 0\f$ and are
/// homogeneous functions of degree, that is
/// \f$\ell\f$ with \f$f(\lambda x, \lambda y, \lambda z) = \lambda^\ell f(x, y, z)\f$
/// for some scaling \f$\lambda\f$.
///
/// TODO make a change of variables u = x + iy, v = x - iy, express angular
/// momentum operators interms of the new variables including L+-. Describe the
/// LRL vector (Laplace-Runge-Lenz vector) in classical orbital mechanics and
/// its application to QM, the extended SO(4) symmetry and raising and lowering
/// operators.
/// Refrences: arxiv.org/pdf/2305.18229

#ifndef HYDROGEN_ATOM_H
#define HYDROGEN_ATOM_H

#include <vector>

#include <glm/glm.hpp>

#include "quantum/qstate.h"

/// Impliments the QState interface for the hydrogen atom.
/// @ingroup quantum_mechanics hydrogen_atom
class HydrogenAtom : public QState {
public:
    struct QuantumNumbers {
        std::size_t m_n;
        std::size_t m_l;
        int m_m;

        QuantumNumbers();
        QuantumNumbers(std::size_t n, std::size_t l, int m);
        bool isValid() const;
        bool operator==(const QuantumNumbers& other) const;
    };

    HydrogenAtom(std::size_t Z, double mu);

    double psi_0_max() const override;
    double psi_n(glm::dvec3 position, std::size_t energy_level_index) const override;
    glm::dvec3 grad_psi_n(glm::dvec3 position, std::size_t energy_level_index) const override;
    std::size_t get_index_from_quantum_numbers(std::vector<int> qn) override;
    std::string get_state_string() const override;
    void find_energy_levels() override;

    double get_energy_eigenvalue(QuantumNumbers quantum_numbers) const;

private:
    std::vector<QuantumNumbers> m_energy_levels_QN;    
    double m_norm;
};

#endif
