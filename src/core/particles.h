/// @file particles.h
/// @brief Classes and structs to setup a particle system.
/// The required code to setup a generic particle system.

#ifndef PARTICLES_H
#define PARTICLES_H

#include <glm/fwd.hpp>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include "core/looplog.h"
#include "core/model.h"

/// Represents a particle in 3D.
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
class Particles {
public:
    /// The particles model.
    Model m_model;
    /// The particles mass.
    double m_mass=1.;

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
    void update(float dt);
    /// Create a new live particle.
    /// If there are avaliable dead particles in the pool, one will be revived.
    /// If there are no avaliable dead particles in the pool, a new live
    /// Particle will be created.
    /// @param decay_rate The particles decay rate.
    /// @param position The position of the new particle.
    void spawn_particle(double decay_rate, glm::vec3 position);
    /// Kill a live particle in the pool.
    /// @param index_in_live The index of the particle in the list of indices of
    /// live particles \ref m_live_IDs.
    void kill_particle(std::size_t index_in_live);
private:
    LoopLog* m_loopLog;
};

#endif
