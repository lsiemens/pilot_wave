#version 330 core

in float view_distance;
in vec3 fragmentColor;

out vec3 color;

// TODO get background color and max view distance
const vec3 background = vec3(0.6, 0.65, 0.7);
const float distance_scale = 20.0;

void main() {
    float weight = exp(-view_distance/distance_scale);
    color = fragmentColor*weight + (1 - weight)*background;
}
