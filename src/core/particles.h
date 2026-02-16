/// @defgroup visualize_distributions Visualize Distributions
/// @brief How to use particle systems to visualize probability distributions.
///
/// The goal here is to visualize a probability distribution with a particle
/// system. This can be done directly by sampling from the probability
/// distribution or some what indirectly by using particles with a variable
/// decay rate. Both methods should work, but I expect that the second method
/// has advantages when dealing with time dependent distributions and
/// distributions with a probability current.
///
/// @section math1 Particle Decay
/// Each particle has a decay rate \f$\Gamma\f$ in Hertz, and label the time
/// step between frames as \f$\Delta t\f$ in seconds. First lets consider the
/// practical implimentation of this decay rate.
///
/// For a particle at the start of the update we compute the probability that
/// the particle decays \f$p = \Gamma \Delta t\f$. Note that while we call this a
/// probability strictly speaking it is not normalized to the interval
/// \f$[0, 1]\f$ as is required for a true probability, this will be addressed
/// shortly. A unit uniform random number is generated and if the random number
/// is less than \f$p\f$ the particle is killed, otherwise it lives till the
/// next step. So the probability that the particle survives till the next step
/// is \f$1 - \min(\Gamma\Delta t, 1)\f$. To simplify this lets assume that the
/// \f$p \ll 1\f$, then the probability of surviving is \f$1 - \Gamma\Delta
/// t\f$. Now lets consider the probabilit that the particle survives for
/// \f$t\f$ secconds, where \f$t = n\Delta t\f$, this probability is
/// \f$P = (1 - \Gamma\Delta t)^n\f$. Note that \f$\Delta t = \frac{t}{n}\f$, so
/// \f$P = \left(1 + \frac{-\Gamma t}{n}\right)^n\f$. Now in the limit that
/// \f$n\f$ is large we find that
/// \f$P = \lim_{n \to \infty} \left(1 + \frac{-\Gamma t}{n}\right)^n = e^{-\Gamma t}\f$.
/// So as implimented, the particles approximate exponential decay based on the
/// provide decay rate.
///
/// So if probability of decaying in a given step is small then this models
/// exponential decay, and the justification for treating the unnormalized
/// quantity \f$\Gamma \Delta t\f$ as a probability is that, given exponential
/// decay the probaility of decaying during a time step is
/// \f$1 - e^{-\Gamma\Delta t} = 1 - (1 - \Gamma\Delta t + \cdots) \approx \Gamma\Delta t\f$.
/// The probability \f$p\f$ is the first nonzero term in the taylor expansion of
/// the exact probability of decaying during the step. When this probability is
/// not small there will be some missmatch between the true value and the
/// approximation but this will only be relavent for particles that have decay
/// rates on the order of the frame rate (life times on order of the duration
/// between frames). That these particles decay marginally faster will likely
/// have almost zero noticable effect.
///
/// So given a decay rate of \f$\Gamma\f$ the probability of the particle
/// surviving \f$t\f$ seconds is \f$P(t) = e^{-\Gamma t}\f$. The expected
/// lifetime of the particle is \f$\bar{t} = E[t] = \frac{1}{\Gamma}\f$. Lastly
/// the half-life of the particle is \f$t_{1/2} = \ln(2) \bar{t}\f$.
///
/// @section math2 Uniform Particle Density
/// In a given volume of space define the number density of particles as \f$n\f$
/// and assume that all of the particles have an identical decay rate of
/// \f$\Gamma\f$, additionaly let new particle be created uniformaly at a rate
/// \f$\lambda\f$ particles per second per unit volume. Then the number density
/// of particles in the volume is described by the differential equation,
/// \f$\frac{d n(t)}{dt} = -\Gamma n(t) + \lambda\f$.
///
/// This differentail equation has the general solution
/// \f$n(t) = Ae^{-\Gamma t} + B\f$. Lets define the initial particle density as
/// \f$n_0 = n(0)\f$ and the density \f$n_\infty = \lim_{t \to \infty} n(t)\f$.
/// Then \f$n_\infty = B\f$ and \f$0 = -\Gamma n_\infty + \lambda\f$ so
/// \f$n_\infty = \frac{\lambda}{\Gamma}\f$. Alternativly at \f$t = 0\f$ we have
/// \f$-\Gamma n_0 + \lambda = -\Gamma A\f$ so \f$A = n_0 - \frac{lambda}{\Gamma} = n_0 - n_\infty\f$.
/// So the general solution if \f$n(t) = (n_0 - n_\infty)e^{-\Gamma t} + n_\infty\f$.
///
/// From the general solution of \f$n(t)\f$, we can see that the number density
/// will asyptocically approch the value \f$n_\infty = \frac{\lambda}{\Gamma}\f$
/// at an exponential rate with half-life \f$t_{1/2} = \frac{\ln(2)}{\Gamma}\f$.
/// 
/// @section math3 Visualizing Probability Densities
/// One way to visualiz a probability density is to generate particles such that
/// the local number density of the particles is proportional to the probability
/// density. Provided the particles do not decay, have a fixed life time, or
/// have a uniform decay rate then the local number density can be set by
/// sampling space weighted by the target probability density, replacing
/// particles as the die from the same distribution.
///
/// Alternatifly given uniform spatial sampling, a variable decay rate can also
/// be used to approach a given target distribution. Start with the following
/// parameters, the maximum mean lifetime \f$\bar{t}_\max\f$, the target total
/// number of particles \f$N\f$ and a normalized probability density \f$\rho(x)\f$.
/// From these parametes we can calculate a uniform and constant particle
/// creation rate \f$\lambda\f$ and a local decay rate \f$\Gamma(x)\f$ such that
/// the asymptotic particle density is proportional to the probability density
/// \f$\rho(x)\f$.
///
/// Define the maximum of the probability density as \f$\rho_\max\f$ and note
/// that the target particle density is \f$n(x) = N\rho(x)\f$. Since the
/// particle density generation rate is constant we can directly calculate it at
/// the point where the desisty is maximized using parameters on hand. Combining
/// the equations for the expected mean lifetime and asymptotic particle density
/// the particle creation rate can be found. Solving the equation gives
/// \f$\lambda = \frac{N\rho_\max}{\bar{t}_\max}\f$, note as expected this is a particle
/// density rate. Now we can use the particle generation rate in the equation
/// for the particle decay rate to find the decay rate as a function of
/// position. After some cansellation
/// \f$\Gamma(x) = \frac{\rho_\max}{\rho \bar{t}_\max}\f$. This decay rate and
/// particle creation rate combine to reproduce the desired target asyptotic
/// distribution. Starting from some arbitrary initial particle density
/// distribution (often it will be n_0(x) = 0) the distribution will approch
/// \f$n(x) = N\rho(x)\f$ at an exponential rate with a maximum half-life of
/// \f$\bar{t}_\max \ln(2)\f$.
///
/// @section math4 Final Notes
/// As noted previousely I expect that using a variable decay rate will have
/// advantages over directly sampliling the probability distribution. When the
/// probability current is zero the distibution is static. Under these
/// contitions all points remain in place and so for every particle the
/// probability density at their location remains fixed. If the probability
/// current is nonzero, then the particles should flow along with the current
/// perserving the expected time dependent probability distrbution from the wave
/// function. Note that for the decaying case, the decay rate of the particles
/// is not static but must continually adjust to match the equations found
/// earlier. Assuming that the particle trajectories could be found exactly,
/// then the direct sampling method is much simpler. However in practice there
/// will be numerical error during the integration of the trajectories and this
/// will lead to a numerical diffusion causing the particles to drift from their
/// expected paths. In the case of direct sampling there is no mechanism to
/// activly correct for these accumulating errors and so over time the
/// distribution can diverge from the true solution. In the case of the decaying
/// test particles, when a given particle moves from the ideal solution due to
/// error it will readjust its properties as if it was generated at that point
/// and since the diribution decays towards the true solution this will add a
/// corrective pressure so whille there may be persistant low level error, it
/// should not be able to compound over time.


#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <functional>

#include <glm/glm.hpp>

#include "core/looplog.h"
#include "core/model.h"
#include "core/rng.h"
#include "core/math_util.h"

/// Represents a particle in 3D.
/// @ingroup visualize_distributions
struct Particle {
    /// The particles ID should be its index in the array of particles.
    const std::size_t m_ID;

    /// The decay rate for this particle.
    double m_decay_rate = 0.0;

    /// The current age of the particle.
    double m_age = 0.0;

    /// The current position of the particle.
    glm::vec3 m_position = glm::vec3(0.0f);

    Particle(std::size_t ID) : m_ID(ID) {}
    Particle(std::size_t ID, double decay_rate, glm::vec3 position) : m_ID(ID), m_decay_rate(decay_rate), m_age(0.0), m_position(position) {}
};

/// A self contained particle system.
/// @ingroup visualize_distributions
class Particles {
public:
    /// The particles model.
    Model m_model;

    /// A pool of preallocated particles.
    std::vector<Particle> m_particle_set;

    /// The indices of live particles in the pool.
    std::vector<std::size_t> m_live_IDs;

    /// The indices of dead particles in the pool.
    std::vector<std::size_t> m_dead_IDs;

    Particles(Model model, std::size_t pool_size);


    /// Draw all live particles in the pool.
    void drawParticles();

    /// Run an update for all live particles in the pool.
    /// @param dt size of the current time step.
    void update(double dt);

    /// Create a new live particle.
    /// If there are available dead particles in the pool, one will be revived.
    /// If there are no available dead particles in the pool, a new live
    /// instance of Particle will be created and added to the pool.
    /// @param decay_rate The particles decay rate.
    /// @param position The position of the new particle.
    void spawn_particle(double decay_rate, glm::vec3 position);

    /// Kill a live particle in the pool.
    /// @param index_in_live The index of the particle in the list of indices of
    /// live particles m_live_IDs.
    void kill_particle(std::size_t index_in_live);
    
    /// calculate the particle generation rate.
    /// Find the right particle generation rate so the longest lived particles
    /// match the target mean lifetime.
    /// @param t_max The mean lifetime of the longest lived particles in seconds.
    /// @param rho_max The maximum probability density of the distribution in
    /// particles per unit volume.
    /// @param N The target total number of particles when converged.
    /// @returns The generation rate \f$\lambda\f$ in particles per second per
    /// unit volume.
    double generation_rate(double t_max, double rho_max, int N) {
        return static_cast<double>(N)*rho_max/t_max;
    }

    /// Calculate the particle decay rate.
    /// Find the particle decay rate at a point that will produce the target
    /// asymptotic particle density.
    /// @param t_max The mean lifetime of the longest lived particles in seconds.
    /// @param rho_x The probability density at a point in particles per unit volume.
    /// @param rho_max The maximum probability density of the distribution in
    /// particles per unit volume.
    /// @returns The decay rate \f$\Gamma(x)\f$ in Hertz.
    double decay_rate(double t_max, double rho_x, double rho_max) {
        return rho_max/(t_max*rho_x);
    }

private:
    LoopLog* m_loopLog;
    RNG m_rng;
};

#endif
