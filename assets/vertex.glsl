#version 330 core

layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec3 vertexColor;

uniform mat4 ProjectionTransform;
uniform mat4 ViewTransform;
uniform mat4 ModelTransform;
uniform float Age;

float Age_Scale = 1.f;
const float PI = 3.141592;

out vec3 fragmentColor;
out float view_distance;

vec3 palette(float arg) {
    return vec3(smoothstep(0.0f, 0.6f, arg),
                smoothstep(0.2f, 0.8f, arg),
                smoothstep(0.5f, 1.0f, arg));
}

void main() {
    vec4 vertexPosition_viewspace = ViewTransform*ModelTransform*vec4(vertexPosition_modelspace, 1);

    if (Age < 0) {
        fragmentColor = vertexColor;
    } else {
        float arg = (1 - exp(-0.2*Age));
        float t2 = (1 - exp(-Age));

        fragmentColor = pow(t2*palette(arg), vec3(0.85f));
    }

    view_distance = length(vertexPosition_viewspace.xyz);
    gl_Position = ProjectionTransform*vertexPosition_viewspace;
}
