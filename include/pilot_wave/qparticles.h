#ifndef QPARTICLES_H
#define QPARTICLES_H

#include <memory>

#define GLEW_STATIC
#include <GL/glew.h>

#include <glm/glm.hpp>

#include "core/looplog.h"
#include "core/particles.h"
#include "quantum/qstate.h"


class QParticles {
public:
    static constexpr float m_particle_size = 0.01f; /// as a fraction of the width
    static constexpr std::size_t m_target_particle_num = 2000;
    static constexpr double m_max_lifetime = 10.;

    std::unique_ptr<QState> m_qstate_uptr;
    std::unique_ptr<Particles> m_particles_uptr;
    double m_rho_0_max;
    double m_spawn_rate;

    QParticles(std::unique_ptr<QState> qstate, GLuint shaderID);
    
    void update(double dt);
    void draw();
private:
    RNG m_rng;
    LoopLog* m_loopLog;
};

#endif
