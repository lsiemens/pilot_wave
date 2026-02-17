#ifndef CONTROLS_H
#define CONTROLS_H

#include <memory>
#include <array>
#include <string>

// keep this before all other OpenGL libraries
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "core/looplog.h"
#include "core/camera.h"
#include "pilot_wave/qparticles.h"

class Controls {
public:
    Controls(GLFWwindow* window, std::shared_ptr<Camera> camera_sprt, std::shared_ptr<QParticles> qparticles_sptr);

    void update(double dt);
    void command(std::string);

private:
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void char_callback(GLFWwindow* window, unsigned int codepoint);
    static void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos);

    std::array<bool, GLFW_KEY_LAST + 1> m_key_pressed = {false};
    std::array<bool, GLFW_KEY_LAST + 1> m_key_held = {false};
    glm::dvec2 m_cursor_pos;
    std::string m_text;
    bool m_is_text_mode = false;

    std::shared_ptr<Camera> m_camera_sptr;
    std::shared_ptr<QParticles> m_qparticles_sptr;
    LoopLog* m_loopLog;
};

#endif
