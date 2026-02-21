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
    static constexpr double m_speed = 3.;
    static constexpr double m_mouseSensitivity = 0.001;


    Controls(GLFWwindow* window, std::shared_ptr<Camera> camera_sprt, std::shared_ptr<QParticles> qparticles_sptr);

    void update(double dt);
    void command(std::string command_str);

private:
    static void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos);
    static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void char_callback(GLFWwindow* window, unsigned int codepoint);

    std::array<bool, GLFW_MOUSE_BUTTON_LAST + 1> m_mouse_button_held = {false};
    std::array<bool, GLFW_KEY_LAST + 1> m_key_pressed = {false};
    std::array<bool, GLFW_KEY_LAST + 1> m_key_held = {false};
    glm::dvec2 m_cursor_pos;
    glm::dvec2 m_cursor_pos_previous;
    glm::dvec2 m_view_angle;
    std::string m_text;
    bool m_is_text_mode = false;
    GLFWwindow* m_window;
    // This variable is reserved for letting the mouse position callback know if
    // it is the first position callback since left clicking. The mouse is
    // enabled/disabled based on the left mouse button, this causes the mouse
    // position to update on release (jump to where it was before clicking). The
    // mouse position call back needs to know that you just left clicked so it
    // can reset m_cursor_pos_previous to be the current position so that
    // differences in the current vs previous pos is continuous.
    bool m_pressed_left_button;

    std::shared_ptr<Camera> m_camera_sptr;
    std::shared_ptr<QParticles> m_qparticles_sptr;
    LoopLog* m_loopLog;
};

#endif
