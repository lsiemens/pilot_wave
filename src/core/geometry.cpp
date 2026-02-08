#include "geometry.h"

#include <cmath>
#include <vector>
#include <stdexcept>

#include <glm/glm.hpp>

constexpr double pi = 3.14159265358979323846;

/// Get the 2d coordinate of the nth vertex in the triangulation of a grid
glm::vec2 meshgrid(int x_res, int y_res, int vertex_index) {
    int triangle_index = vertex_index / 3; // index of the triangle
    int triangle_vertex_index = vertex_index % 3; // label of the vertex within the triangle

    int quad_index = triangle_index / 2; // index of the quad
    int quad_triangle_index = triangle_index % 2; // label of the triangle within the quad

    int column_index = quad_index / (x_res - 1); // index of the column
    int row_index = quad_index % (x_res - 1); // index of the row

    // lower left coordinate of triangle on the unit square
    glm::vec2 unit_pos;
    unit_pos.x = static_cast<float>(row_index)/static_cast<float>(x_res - 1);
    unit_pos.y = static_cast<float>(column_index)/static_cast<float>(y_res - 1);

    // offset for the upper right vertex
    if (triangle_vertex_index == 1) {
        unit_pos.x += 1.0f/static_cast<float>(x_res - 1);
        unit_pos.y += 1.0f/static_cast<float>(y_res - 1);
    }

    // offsets for the off diagonal vertex. The spesific offset
    // is dependent on if it is the upper or lower triangle
    if (triangle_vertex_index == 2){
        if (quad_triangle_index == 0) {
            unit_pos.x += 1.0f/static_cast<float>(x_res - 1);
        } else {
            unit_pos.y += 1.0f/static_cast<float>(y_res - 1);
        }
    }
    return unit_pos;
}

Model make_model(SurfaceMap surface_map, SurfaceMap color_map, int x_res, int y_res, GLuint shaderID) {
    int num_vertices = (x_res - 1)*(y_res - 1)*2*3; // (x_res - 1)*(y_res - 1) quads, 2 triangles per quad, 3 points per triangle

    std::vector<GLfloat> g_vertex_buffer(num_vertices*3);
    std::vector<GLfloat> g_color_buffer(num_vertices*3);

    for (int vertex_index = 0; vertex_index < num_vertices; vertex_index++) {
        glm::vec2 unit_pos = meshgrid(x_res, y_res, vertex_index);

        glm::vec3 vertex_pos = surface_map(unit_pos);
        glm::vec3 vertex_col = color_map(unit_pos);

        g_vertex_buffer[0 + 3*vertex_index] = vertex_pos.x;
        g_vertex_buffer[2 + 3*vertex_index] = vertex_pos.y;
        g_vertex_buffer[1 + 3*vertex_index] = vertex_pos.z;

        g_color_buffer[0 + 3*vertex_index] = vertex_col.r;
        g_color_buffer[1 + 3*vertex_index] = vertex_col.g;
        g_color_buffer[2 + 3*vertex_index] = vertex_col.b;
    }

    Model model = Model(shaderID);
    model.setVertexBuffer(g_vertex_buffer.data(), static_cast<GLsizei>(g_vertex_buffer.size())*sizeof(GLfloat));
    model.setColorBuffer(g_color_buffer.data(), static_cast<GLsizei>(g_color_buffer.size())*sizeof(GLfloat));
    return model;
}

/// Insert three glm::vec3 defining the values at the vertices of a triangle into a buffer.
void insert_triangle(std::vector<GLfloat>& buffer, std::size_t offset, glm::vec3 vert_1, glm::vec3 vert_2, glm::vec3 vert_3) {
    if (buffer.size() < 9 + 3*offset) {
        throw std::runtime_error("Inserting triangle into buffer. Index out of bounds!");
    }

    buffer[0 + 3*offset] = vert_1.x;
    buffer[2 + 3*offset] = vert_1.y;
    buffer[1 + 3*offset] = vert_1.z;

    buffer[3 + 3*offset] = vert_2.x;
    buffer[5 + 3*offset] = vert_2.y;
    buffer[4 + 3*offset] = vert_2.z;

    buffer[6 + 3*offset] = vert_3.x;
    buffer[8 + 3*offset] = vert_3.y;
    buffer[7 + 3*offset] = vert_3.z;
}

Model make_tetrahedron(float size, GLuint shaderID) {
    int num_vertices = 4*3; // (x_res - 1)*(y_res - 1) quads, 2 triangles per quad, 3 points per triangle

    size *= std::sqrt(8.f/3.f);

    std::vector<GLfloat> g_vertex_buffer(num_vertices*3);
    std::vector<GLfloat> g_color_buffer(num_vertices*3);

    glm::vec3 vert_A = size*glm::vec3( 0.5f,-1/(2*std::sqrt(3.f)),-1/(2*std::sqrt(6.f)));
    glm::vec3 vert_B = size*glm::vec3( 0.0f, 1/std::sqrt(3.f),    -1/(2*std::sqrt(6.f)));
    glm::vec3 vert_C = size*glm::vec3(-0.5f,-1/(2*std::sqrt(3.f)),-1/(2*std::sqrt(6.f)));
    glm::vec3 vert_D = size*glm::vec3( 0.0f, 0.0f,                 std::sqrt(3.f/8.f));

    glm::vec3 col_A = glm::vec3( 1.f, 0.f, 0.f);
    glm::vec3 col_B = glm::vec3( 0.f, 0.f, 1.f);
    glm::vec3 col_C = glm::vec3( 0.f, 1.f, 0.f);
    glm::vec3 col_D = glm::vec3( 1.f, 1.f, 1.f);

    insert_triangle(g_vertex_buffer, 0, vert_A, vert_C, vert_B);
    insert_triangle(g_vertex_buffer, 3, vert_A, vert_B, vert_D);
    insert_triangle(g_vertex_buffer, 6, vert_B, vert_C, vert_D);
    insert_triangle(g_vertex_buffer, 9, vert_C, vert_A, vert_D);

    insert_triangle(g_color_buffer, 0, col_A, col_C, col_B);
    insert_triangle(g_color_buffer, 3, col_A, col_B, col_D);
    insert_triangle(g_color_buffer, 6, col_B, col_C, col_D);
    insert_triangle(g_color_buffer, 9, col_C, col_A, col_D);

    Model model = Model(shaderID);
    model.setVertexBuffer(g_vertex_buffer.data(), static_cast<GLsizei>(g_vertex_buffer.size())*sizeof(GLfloat));
    model.setColorBuffer(g_color_buffer.data(), static_cast<GLsizei>(g_color_buffer.size())*sizeof(GLfloat));
    return model;
}

/// The embedding of a sphere into R^3.
glm::vec3 sphere_map(float radius, glm::vec2 uv_position) {
    glm::vec3 position;

    float theta = 2*static_cast<float>(pi)*uv_position.x;
    float phi = static_cast<float>(pi)*uv_position.y;

    position.x = radius*std::sin(phi)*std::cos(theta);
    position.y = radius*std::sin(phi)*std::sin(theta);
    position.z = radius*std::cos(phi);

    return position;
}

/// The embedding of a torus into R^3.
glm::vec3 torus_map(float minor_radius, float major_radius, glm::vec2 uv_position) {
    glm::vec3 position;

    float theta = 2*static_cast<float>(pi)*uv_position.x;
    float phi = 2*static_cast<float>(pi)*uv_position.y;

    position.x = (major_radius + minor_radius*std::cos(theta))*std::cos(phi);
    position.y = (major_radius + minor_radius*std::cos(theta))*std::sin(phi);
    position.z = minor_radius*std::sin(theta);

    return position;
}

/// A default uv color map.
glm::vec3 color_map(glm::vec2 uv_position) {
    glm::vec3 color;
    color.r = uv_position.x;
    color.g = uv_position.y;
    color.b = 0.f;
    return color;
}

Model make_sphere(float radius, int x_res, int y_res, GLuint shaderID) {
    return make_model([radius](glm::vec2 uv) -> glm::vec3 {return sphere_map(radius, uv);}, color_map, x_res, y_res, shaderID);
}

Model make_torus(float minor_radius, float major_radius, int x_res, int y_res, GLuint shaderID) {
    return make_model([minor_radius, major_radius](glm::vec2 uv) -> glm::vec3 {return torus_map(minor_radius, major_radius, uv);}, color_map, x_res, y_res, shaderID);
}
