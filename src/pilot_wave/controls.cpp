#include "controls.h"

Controls::Controls(GLFWwindow* window, std::shared_ptr<Camera> camera_sptr, std::shared_ptr<QParticles> qparticles_sptr) {
    assert(camera_sptr);
    assert(qparticles_sptr);
    m_camera_sptr = camera_sptr;
    m_qparticles_sptr = qparticles_sptr;
    m_loopLog = LoopLog::getInstance();

    // register callbacks
    glfwSetWindowUserPointer(window, this);
    glfwSetKeyCallback(window, key_callback);
    glfwSetCharCallback(window, char_callback);
    glfwSetCursorPosCallback(window, cursor_pos_callback);
}

void Controls::update(double dt) {
    double horizontalAngle = 3.13, verticalAngle = 0.0;
    float speed = 3.f, mouseSensitivity = 0.001f;

    horizontalAngle = -mouseSensitivity*m_cursor_pos.x;
    verticalAngle = mouseSensitivity*m_cursor_pos.y;

    glm::vec3 direction = glm::vec3(std::cos(verticalAngle)*std::sin(horizontalAngle),
                     std::sin(verticalAngle),
                     std::cos(verticalAngle)*std::cos(horizontalAngle));

    glm::vec3 right = glm::vec3(std::sin(horizontalAngle - 3.14f/2.f),
                     0.f,
                     std::cos(horizontalAngle - 3.14/2.f));

    glm::vec3 up = glm::cross(right, direction);
    glm::vec3 delta_position = glm::vec3(0.f, 0.f, 0.f);

    if (m_key_held[GLFW_KEY_W]) {
        delta_position = static_cast<float>(dt)*direction*speed;
    }

    if (m_key_held[GLFW_KEY_S]) {
        delta_position = -static_cast<float>(dt)*direction*speed;
    }

    if (m_key_held[GLFW_KEY_D]) {
        delta_position = static_cast<float>(dt)*right*speed;
    }

    if (m_key_held[GLFW_KEY_A]) {
        delta_position = -static_cast<float>(dt)*right*speed;
    }

    if (m_key_held[GLFW_KEY_R]) {
        delta_position = static_cast<float>(dt)*up*speed;
    }

    if (m_key_held[GLFW_KEY_F]) {
        delta_position = -static_cast<float>(dt)*up*speed;
    }

    m_camera_sptr->m_direction = direction;
    m_camera_sptr->m_up = up;
    m_camera_sptr->m_position += delta_position;

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
}

void Controls::command(std::string) {
    //Process command
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

    if ((key > 0) and (key < GLFW_KEY_LAST + 1)) {
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
}
