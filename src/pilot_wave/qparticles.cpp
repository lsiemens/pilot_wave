#include "pilot_wave/qparticles.h"

#include <iostream>

#include "core/geometry.h"
#include "core/model.h"
#include "core/math_util.h"

QParticles::QParticles(std::unique_ptr<QState> qstate, GLuint shaderID)
        : m_qstate_uptr(std::move(qstate)) {

    double width = m_qstate_uptr->get_width();
    Model model = make_tetrahedron(m_particle_size*width, shaderID);

    m_particles_uptr = std::make_unique<Particles>(model, m_target_particle_num);

    m_rho_0_max = m_qstate_uptr->probability_density_0_max();

    double volume = width*width*width;
    m_spawn_rate = m_particles_uptr->generation_rate(m_max_lifetime, m_rho_0_max, m_target_particle_num)*volume;

    m_loopLog = LoopLog::getInstance();
    m_rng = RNG();
}

void QParticles::update(double dt) {
    double width = m_qstate_uptr->get_width();
    glm::vec3 origin = m_qstate_uptr->get_origin();

    double lambda = dt*m_spawn_rate;
    std::size_t num_new_particles = m_rng.poisson(lambda);
    for (std::size_t i = 0; i < num_new_particles; i++) {
        glm::vec3 offset;
        offset.x = width*(m_rng.uniform() - 0.5);
        offset.y = width*(m_rng.uniform() - 0.5);
        offset.z = width*(m_rng.uniform() - 0.5);
        m_particles_uptr->spawn_particle(0.f, origin + offset);
    }
    m_loopLog->m_log << "QParticles | dt*spawn_rate: [" << lambda << "]\n";
    m_loopLog->m_log << "QParticles | new particles: [" << num_new_particles << "]\n";

    m_qstate_uptr->update(dt);
    m_particles_uptr->update(dt);

    float dt_f = static_cast<float>(dt);

    for (std::size_t i = m_particles_uptr->m_live_IDs.size(); i-- > 0; ) {
        Particle& particle = m_particles_uptr->m_particle_set[m_particles_uptr->m_live_IDs[i]];
        particle.m_position += dt_f*glm::vec3(0.f, 0.f, 0.f);

        double rho = m_qstate_uptr->probability_density(particle.m_position);
        glm::dvec3 v_rho = m_qstate_uptr->probability_current(particle.m_position);

        particle.m_decay_rate = m_particles_uptr->decay_rate(m_max_lifetime, rho, m_rho_0_max);
    }
}

void QParticles::draw() {
    m_particles_uptr->drawParticles();
}
