#include "core/particles.h"

#include <cassert>
#include <iostream>
#include <algorithm>

// /todo remove
#include <cstdlib>
#include <ctime>

Particles::Particles(Model model, std::size_t pool_size) : m_model(model) {
    m_particle_set.reserve(pool_size);
    m_dead_IDs.reserve(pool_size);


    // /todo remove
    std::srand(std::time(0));

    for (std::size_t i = 0; i < pool_size; i++) {
        m_particle_set.emplace_back(i);
        m_dead_IDs.emplace_back(i);
    }

    m_loopLog = LoopLog::getInstance();
}

void Particles::drawParticles() {
    for (std::size_t i = 0; i < m_live_IDs.size(); i++) {
        std::size_t particle_ID = m_live_IDs[i];

        // draw particle given by particle_ID
        this->m_model.drawModel(glm::translate(glm::mat4(1.f), m_particle_set[particle_ID].m_position));
    }
}

void Particles::update(float dt) {
    m_loopLog->m_log << "Particle system | pool size:" << m_particle_set.size() << "\n";
    m_loopLog->m_log << "\tparticles state [live | dead]: [" << m_live_IDs.size() << " | " << m_dead_IDs.size() << "]\n";

    /// When updating iterate over live_IDs in reverse, this way objects can be
    /// removed by swapping with the last element and then popping with out
    /// effecting what element should be processed next.
    for (std::size_t i = m_live_IDs.size(); i-- > 0; ) {
        std::size_t particle_ID = m_live_IDs[i];

        // update particle given by particle_ID
        m_particle_set[particle_ID].m_age += dt;
        m_particle_set[particle_ID].m_position += glm::vec3(0.0f, -dt*1.0f, 0.0f);

        if (m_particle_set[particle_ID].m_age > 100.0) {
            kill_particle(i);
            continue;
        } 
    }
}

void Particles::spawn_particle(double decay_rate, glm::vec3 position) {
    std::size_t particle_ID;

    // Use dead particles if they are avaliable, else create new particles
    if (m_dead_IDs.size() > 0) {
        particle_ID = m_dead_IDs.back();
        m_particle_set[particle_ID].m_decay_rate = decay_rate;
        m_particle_set[particle_ID].m_age = 0.0;
        m_particle_set[particle_ID].m_position = position;


        // Revive particle
        m_dead_IDs.pop_back();
        m_live_IDs.emplace_back(particle_ID);
    } else {
        assert(m_live_IDs.size() == m_particle_set.size());
        particle_ID = m_live_IDs.size();

        // Create particle
        m_particle_set.emplace_back(particle_ID, decay_rate, position);
        m_live_IDs.emplace_back(particle_ID);
    }

    // /todo remove
    float max = RAND_MAX;
    float x = 2*(std::rand()/max - 0.5);
    float y = 2*(std::rand()/max - 0.5);
    m_particle_set[particle_ID].m_position = glm::vec3(100*x, 10.0f, 100*y);
}

void Particles::kill_particle(std::size_t index_in_live) {
    std::size_t particle_ID = m_live_IDs[index_in_live];

    std::swap(m_live_IDs[index_in_live], m_live_IDs.back());
    assert(m_live_IDs.back() == particle_ID);

    // Kill particle
    m_live_IDs.pop_back();
    m_dead_IDs.emplace_back(particle_ID);
}
