#ifndef PHYSICS_UTIL_H
#define PHYSICS_UTIL_H

#include <functional>

#include <glm/glm.hpp>

/// Vector field type
/// Defines a vector filed as a std::function.
/// @param position The position of a point.
/// @param t_offset An offset from the current time in secconds.
/// @returns The velocity at the point as a glm::vec3.
using VectorField = std::function<glm::vec3(glm::vec3 position, double t_offset)>;

#endif
