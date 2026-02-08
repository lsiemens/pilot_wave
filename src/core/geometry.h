#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <functional>

#define GLEW_STATIC
#include <GL/glew.h>

#include "core/model.h"

/// A embedding of a surface.
/// @param uv_position A glm::vec2 vector of uv coordinates.
/// @returns A glm::vec3 position vector of a point on the surface.
using SurfaceMap = std::function<glm::vec3(glm::vec2 uv_position)>;

/// Generate a model from surface maps.
/// @param surface_map A map from uv coordinates to the surface.
/// @param color_map A map from uv coordinates to vertex colors.
/// @param x_res The u resolution of the surface.
/// @param y_res The v resolution of the surface.
/// @param shaderID The id of a shader to apply;
Model make_model(SurfaceMap surface_map, SurfaceMap color_map, int x_res, int y_res, GLuint shaderID);

/// Generate a tetrahedron.
/// @param size The radial distance from the center to all vertices.
/// @param shaderID The id of a shader to apply;
Model make_tetrahedron(float size, GLuint shaderID);

/// Generate a sphere.
/// @param radius The radius of the sphere.
/// @param x_res The longditudinal resolution.
/// @param y_res The latitudinal resolution.
/// @param shaderID The id of a shader to apply;
Model make_sphere(float radius, int x_res, int y_res, GLuint shaderID);

/// Generate a torus.
/// @param minor_radius The radius of the minor direction.
/// @param major_radius The radius of the major direction.
/// @param x_res The resolution about the minor direction.
/// @param y_res The resolution about the major direction.
/// @param shaderID The id of a shader to apply;
Model make_torus(float minor_radius, float major_radius, int x_res, int y_res, GLuint shaderID);

#endif
