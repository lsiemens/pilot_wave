#include "pilot_wave/controls.h"

Controls::Controls(GLFWwindow* window, std::shared_ptr<Camera> camera_sptr,
                   std::shared_ptr<QParticles> qparticles_sptr) : m_window(window) {
    assert(camera_sptr);
    assert(qparticles_sptr);
    m_camera_sptr = camera_sptr;
    m_qparticles_sptr = qparticles_sptr;
    m_loopLog = LoopLog::getInstance();

    // Input settings
    glfwSetInputMode(m_window, GLFW_STICKY_KEYS, GL_TRUE);
    glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    // register callbacks
    glfwSetWindowUserPointer(m_window, this);
    glfwSetMouseButtonCallback(m_window, mouse_button_callback);
    glfwSetKeyCallback(m_window, key_callback);
    glfwSetCharCallback(m_window, char_callback);
    glfwSetCursorPosCallback(m_window, cursor_pos_callback);
}

void Controls::update(double dt) {
    if (m_mouse_button_held[GLFW_MOUSE_BUTTON_LEFT]) {
        m_view_angle += m_mouseSensitivity*(m_cursor_pos - m_cursor_pos_previous);
    }

    glm::dvec3 direction = glm::vec3(std::cos(m_view_angle.y)*std::sin(m_view_angle.x),
                                    std::sin(m_view_angle.y),
                                    std::cos(m_view_angle.y)*std::cos(m_view_angle.x));

    glm::dvec3 right = glm::vec3(std::sin(m_view_angle.x - 3.14/2.),
                                0.,
                                std::cos(m_view_angle.x - 3.14/2.));

    glm::dvec3 up = glm::cross(right, direction);
    glm::dvec3 delta_position = glm::vec3(0., 0., 0.);

    if (m_key_held[GLFW_KEY_W]) {
        delta_position += direction;
    }

    if (m_key_held[GLFW_KEY_S]) {
        delta_position += -direction;
    }

    if (m_key_held[GLFW_KEY_D]) {
        delta_position += right;
    }

    if (m_key_held[GLFW_KEY_A]) {
        delta_position += -right;
    }

    if (m_key_held[GLFW_KEY_R]) {
        delta_position += up;
    }

    if (m_key_held[GLFW_KEY_F]) {
        delta_position += -up;
    }

    if (glm::length(delta_position) > 0.) {
        delta_position = glm::normalize(delta_position);
    }

    m_camera_sptr->m_direction = direction;
    m_camera_sptr->m_up = up;
    m_camera_sptr->m_position += dt*m_speed*delta_position;

    // Ingest key press, it will be reset after this point
    if (m_key_pressed[GLFW_KEY_E]) {
        std::size_t energy_level = m_qparticles_sptr->m_qstate_uptr->get_energy_level();
        m_qparticles_sptr->m_qstate_uptr->set_energy_level(energy_level + 1);

        m_key_pressed[GLFW_KEY_E] = false;
    }

    // Ingest key press, it will be reset after this point
    if (m_key_pressed[GLFW_KEY_Q]) {
        std::size_t energy_level = m_qparticles_sptr->m_qstate_uptr->get_energy_level();
        if (energy_level > 0) {
            m_qparticles_sptr->m_qstate_uptr->set_energy_level(energy_level - 1);
        }
        m_key_pressed[GLFW_KEY_Q] = false;
    }

    if (m_is_text_mode) {
        m_loopLog->m_log << "Controls mode [Text | Keys]: [Text]\n";
        m_loopLog->m_log << "\tText: [" << m_text << "]\n";
    } else {
        m_loopLog->m_log << "Controls mode [Text | Keys]: [Keys]\n";
        m_loopLog->m_log << "\tText: [" << m_text << "]\n";
    }
    m_loopLog->m_log << "Mouse input [left | middle | right]: [" << m_mouse_button_held[GLFW_MOUSE_BUTTON_LEFT] << " | " << m_mouse_button_held[GLFW_MOUSE_BUTTON_MIDDLE] << " | " << m_mouse_button_held[GLFW_MOUSE_BUTTON_RIGHT] << "]\n";

    // This should remain at the end of update
    m_cursor_pos_previous = m_cursor_pos;
}

void Controls::command(std::string) {
    //Process command
}

void Controls::mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    Controls* controls_ptr = static_cast<Controls*>(glfwGetWindowUserPointer(window));

    if ((button >= 0) and (button < GLFW_MOUSE_BUTTON_LAST + 1)) {
        if (action == GLFW_PRESS) {
            controls_ptr->m_mouse_button_held[button] = true;
        }

        if (action == GLFW_RELEASE) {
            controls_ptr->m_mouse_button_held[button] = false;
        }
    }

    // capture mouse on left click
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            controls_ptr->m_pressed_left_button = true;
            glfwSetInputMode(controls_ptr->m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        }

        if (action == GLFW_RELEASE) {
            glfwSetInputMode(controls_ptr->m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        }
    }
}

void Controls::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    Controls* controls_ptr = static_cast<Controls*>(glfwGetWindowUserPointer(window));

    // Set control mode
    if ((key == GLFW_KEY_ENTER) and (action == GLFW_PRESS)) {
        if (controls_ptr->m_is_text_mode) {
            controls_ptr->m_is_text_mode = false;
            controls_ptr->command(controls_ptr->m_text);
            controls_ptr->m_text = "";
        } else {
            controls_ptr->m_is_text_mode = true;
        }
    }

    if ((key >= 0) and (key < GLFW_KEY_LAST + 1)) {
        if (action == GLFW_PRESS) {
            // only limit pressing and not release in text mode to keys do not
            // get stuck
            if (not controls_ptr->m_is_text_mode) { 
                controls_ptr->m_key_held[key] = true;
                controls_ptr->m_key_pressed[key] = true;
            }
        }

        if (action == GLFW_RELEASE) {
            controls_ptr->m_key_held[key] = false;
        }
    }

    // Check special keys for text input
    if (controls_ptr->m_is_text_mode) {
        // Only check special characters during text mode
        if ((key == GLFW_KEY_DELETE) or (key == GLFW_KEY_BACKSPACE)) {
            if ((action == GLFW_PRESS) or (action == GLFW_REPEAT)) {
                if (controls_ptr->m_text.size() > 0) {
                    controls_ptr->m_text.pop_back();
                }
            }
        }
    }
}

void Controls::char_callback(GLFWwindow* window, unsigned int codepoint) {
    Controls* controls_ptr = static_cast<Controls*>(glfwGetWindowUserPointer(window));

    if (controls_ptr->m_is_text_mode) {
        controls_ptr->m_text += static_cast<char>(codepoint);
    }
}

void Controls::cursor_pos_callback(GLFWwindow* window, double xpos, double ypos) {
    Controls* controls_ptr = static_cast<Controls*>(glfwGetWindowUserPointer(window));

    controls_ptr->m_cursor_pos.x = xpos;
    controls_ptr->m_cursor_pos.y = ypos;

    if (controls_ptr->m_pressed_left_button) {
        controls_ptr->m_pressed_left_button = false;
        controls_ptr->m_cursor_pos_previous = controls_ptr->m_cursor_pos;
    }
}
