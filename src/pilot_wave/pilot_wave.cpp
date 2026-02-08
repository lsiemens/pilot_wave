#include <cmath>
#include <iostream>

// keep this before all other OpenGL libraries
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "core/looplog.h"
#include "core/frame_timer.h"
#include "core/model.h"
#include "core/object.h"
#include "core/particles.h"
#include "core/camera.h"
#include "core/shaders.h"
#include "core/path_util.h"
#include "core/rng.h"
#include "core/physics_util.h"
#include "core/geometry.h"

void Controlls(float dt, GLFWwindow* window, Camera &camera) {
    double horizontalAngle = 3.13, verticalAngle = 0.0;
    float speed = 3.f, mouseSensitivity = 0.001f;

    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    horizontalAngle = -mouseSensitivity*xpos;
    verticalAngle = mouseSensitivity*ypos;

    glm::vec3 direction = glm::vec3(std::cos(verticalAngle)*std::sin(horizontalAngle),
                     std::sin(verticalAngle),
                     std::cos(verticalAngle)*std::cos(horizontalAngle));

    glm::vec3 right = glm::vec3(std::sin(horizontalAngle - 3.14f/2.f),
                     0.f,
                     std::cos(horizontalAngle - 3.14/2.f));

    glm::vec3 up = glm::cross(right, direction);
    glm::vec3 delta_position = glm::vec3(0.f, 0.f, 0.f);

    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
        delta_position = dt*direction*speed;
    }
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
        delta_position = -dt*direction*speed;
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
        delta_position = dt*right*speed;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) {
        delta_position = -dt*right*speed;
    }
    if (glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS) {
        delta_position = dt*up*speed;
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS) {
        delta_position = -dt*up*speed;
    }

    camera.m_direction = direction;
    camera.m_up = up;
    camera.m_position += delta_position;
}

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

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Initalize shader
    GLuint shaderID = LoadShaders("assets/vertex.glsl", "assets/fragment.glsl");

    AdvancedTimer timer = AdvancedTimer();
    Camera camera = Camera(shaderID);

    Object sphere = Object(make_sphere(1.0f, 100, 100, shaderID));
    sphere.m_position = glm::vec3(3.0f, 0.0f, -3.0f);
    sphere.m_velocity = glm::vec3(0.0f, 10.0f, 0.0f);
    sphere.m_acceleration = glm::vec3(0.0f, -9.81f, 0.0f);

    Object torus = Object(make_torus(0.5f, 1.f, 100, 100, shaderID));
    torus.m_position = glm::vec3(-3.0f, 0.0f, -3.0f);

    RNG rng = RNG();

    VectorField velocity = [](glm::vec3 position, double t_offset) {return glm::vec3(0.f, -1.f, 0.f);};
    Particles test_particles(velocity, make_tetrahedron(.1f, shaderID), 2000);

    float dt;
    do {
        for (int i = 0; i < 300; i++) {
            float x = 100.f*static_cast<float>(2*rng.uniform() - 1.f);
            float y = 100.f*static_cast<float>(2*rng.uniform() - 1.f);
            test_particles.spawn_particle(0.1f, glm::vec3(x, 10.0f, y));
        }

        // Timing
        dt = static_cast<float>(timer.timer());
        loopLog->flush();

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //Camera
        Controlls(dt, window, camera);

        glUseProgram(shaderID);
        camera.update();

        sphere.update(dt);
        torus.update(dt);
        test_particles.update(dt);

        sphere.drawObject();
        torus.drawObject();
        test_particles.drawParticles();

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while ((glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS) && (glfwWindowShouldClose(window) == 0));
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
