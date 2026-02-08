#version 330 core

layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 vertexColor;

uniform mat4 ProjectionTransform;
uniform mat4 ViewTransform;
uniform mat4 ModelTransform;

out vec3 fragmentColor;
out float view_distance;

void main() {
    vec4 vertexPosition_viewspace = ViewTransform*ModelTransform*vec4(vertexPosition_modelspace, 1);

    view_distance = length(vertexPosition_viewspace.xyz);
    fragmentColor = vertexColor;
    gl_Position = ProjectionTransform*vertexPosition_viewspace;
}
