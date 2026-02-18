#include <cmath>
#include <iostream>
#include <memory>

// keep this before all other OpenGL libraries
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "core/looplog.h"
#include "core/frame_timer.h"
#include "core/model.h"
#include "core/object.h"
#include "core/camera.h"
#include "core/shaders.h"
#include "core/path_util.h"
#include "core/geometry.h"

#include "quantum/square_well.h"
#include "pilot_wave/qparticles.h"
#include "pilot_wave/controls.h"

int main() {
    LoopLog* loopLog = LoopLog::getInstance();
    if (!glfwInit()) {
        std::cerr << "Failed to initalize GLFW\n";
        return -1;
    }

    GLFWwindow* window;
    int width=1024, height=768;
    glfwWindowHint(GLFW_SAMPLES, 4); // set multi sampling factor
    window = glfwCreateWindow( width, height, "Pilot Wave", NULL, NULL);
    if (window == NULL) {
        std::cerr << "Failed to create GLFW window.\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(0);
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initalize GLEW.\n";
        return -1;
    }

    glClearColor(.6f, .65f, .7f, 1.f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_MULTISAMPLE); // enable multi sampling

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Initalize shader
    GLuint shaderID = LoadShaders("assets/vertex.glsl", "assets/fragment.glsl");

    AdvancedTimer timer = AdvancedTimer();
    std::shared_ptr<Camera> camera_sptr = std::make_shared<Camera>(shaderID);

    Object sphere = Object(make_sphere(1.0f, 100, 100, shaderID));
    sphere.m_position = glm::vec3(3.0f, 0.0f, -3.0f);
    sphere.m_velocity = glm::vec3(0.0f, 10.0f, 0.0f);
    sphere.m_acceleration = glm::vec3(0.0f, -9.81f, 0.0f);

    Object torus = Object(make_torus(0.5f, 1.f, 100, 100, shaderID));
    torus.m_position = glm::vec3(-3.0f, 0.0f, -3.0f);

    auto qstate_uptr = std::make_unique<SquareWell>(2.);
    std::shared_ptr<QParticles> qparticles_sptr = std::make_shared<QParticles>(std::move(qstate_uptr), shaderID);
    qparticles_sptr->m_qstate_uptr->set_energy_level(1);

    Controls controls = Controls(window, camera_sptr, qparticles_sptr);

    float dt;
    do {
        // Timing
        dt = static_cast<float>(timer.timer());
        loopLog->flush();

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //Camera
        controls.update(dt);

        glUseProgram(shaderID);
        camera_sptr->update();

        sphere.update(dt);
        torus.update(dt);
        qparticles_sptr->update(dt);

        sphere.drawObject();
        torus.drawObject();
        qparticles_sptr->draw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while ((glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS) && (glfwWindowShouldClose(window) == 0));
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
