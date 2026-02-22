#include "core/particles.h"

#include <cassert>
#include <algorithm>

#include <glm/gtx/transform.hpp>

Particles::Particles(Model model, std::size_t pool_size) : m_model(model) {
    m_particle_set.reserve(pool_size);
    m_dead_IDs.reserve(pool_size);

    for (std::size_t i = 0; i < pool_size; i++) {
        m_particle_set.emplace_back(i);
        m_dead_IDs.emplace_back(i);
    }

    m_loopLog = LoopLog::getInstance();
    m_rng = RNG();
}

void Particles::drawParticles() {
    for (std::size_t i = 0; i < m_live_IDs.size(); i++) {
        std::size_t particle_ID = m_live_IDs[i];

        // draw particle given by particle_ID
        glm::mat4 translation = glm::translate(glm::mat4(1.f), m_particle_set[particle_ID].m_position);
        this->m_model.drawModel(translation, static_cast<float>(m_particle_set[particle_ID].m_age));
    }
}

void Particles::update(double dt) {
    m_loopLog->m_log << "Particle system | pool size:" << m_particle_set.size() << "\n";
    m_loopLog->m_log << "\tparticles state [live | dead]: [" << m_live_IDs.size() << " | " << m_dead_IDs.size() << "]\n";

    /// When updating iterate over live_IDs in reverse, this way objects can be
    /// removed by swapping with the last element and then popping with out
    /// effecting what element should be processed next.
    for (std::size_t i = m_live_IDs.size(); i-- > 0; ) {
        // Get current particle
        Particle& particle = m_particle_set[m_live_IDs[i]];

        // Check for decay
        if (particle.m_decay_rate > 0.0) {
            double decay_prob = dt*particle.m_decay_rate;
            double test_value = m_rng.uniform();

            // kill the particle with probability equal to decay_prob
            if (test_value < decay_prob) {
                kill_particle(i);
                continue;
            }
        }

        particle.m_age += dt;
    }
}

void Particles::spawn_particle(double decay_rate, glm::vec3 position) {
    // Use dead particles if they are avaliable, else create new particles
    if (m_dead_IDs.size() > 0) {
        std::size_t particle_ID = m_dead_IDs.back();
        m_particle_set[particle_ID].m_decay_rate = decay_rate;
        m_particle_set[particle_ID].m_age = 0.0;
        m_particle_set[particle_ID].m_position = position;


        // Revive particle
        m_dead_IDs.pop_back();
        m_live_IDs.emplace_back(particle_ID);
    } else {
        assert(m_live_IDs.size() == m_particle_set.size());
        std::size_t particle_ID = m_live_IDs.size();

        // Create particle
        m_particle_set.emplace_back(particle_ID, decay_rate, position);
        m_live_IDs.emplace_back(particle_ID);
    }
}

void Particles::kill_particle(std::size_t index_in_live) {
    std::size_t particle_ID = m_live_IDs[index_in_live];

    std::swap(m_live_IDs[index_in_live], m_live_IDs.back());
    assert(m_live_IDs.back() == particle_ID);

    // Kill particle
    m_live_IDs.pop_back();
    m_dead_IDs.emplace_back(particle_ID);
}
