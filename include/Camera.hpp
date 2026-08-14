#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// Simple orbital camera around a target point.
//
// Controls are driven externally (main.cpp pumps GLFW input into these
// methods). State is kept in spherical coordinates so mouse drags map
// naturally to yaw/pitch rotations.
namespace waveguide {

class OrbitCamera {
public:
    OrbitCamera() = default;

    void setTarget(const glm::vec3& t) { target_ = t; }
    void setDistance(float d)           { distance_ = d; }

    // Mouse drag deltas (in pixels) => rotate around target.
    void orbit(float dxPixels, float dyPixels) {
        const float sensitivity = 0.005f;
        yaw_   += dxPixels * sensitivity;
        pitch_ += dyPixels * sensitivity;

        const float lim = glm::radians(89.0f);
        if (pitch_ >  lim) pitch_ =  lim;
        if (pitch_ < -lim) pitch_ = -lim;
    }

    // Scroll wheel => zoom in/out.
    void zoom(float delta) {
        distance_ *= std::pow(0.9f, delta);
        if (distance_ < 0.005f) distance_ = 0.005f;
        if (distance_ > 100.0f) distance_ = 100.0f;
    }

    glm::vec3 position() const {
        const float cp = std::cos(pitch_);
        const float sp = std::sin(pitch_);
        const float cy = std::cos(yaw_);
        const float sy = std::sin(yaw_);
        return target_ + distance_ * glm::vec3(cp * sy, sp, cp * cy);
    }

    glm::mat4 viewMatrix() const {
        return glm::lookAt(position(), target_, glm::vec3(0.0f, 1.0f, 0.0f));
    }

    glm::mat4 projectionMatrix(float aspect) const {
        return glm::perspective(glm::radians(fovDeg_), aspect, 0.001f, 1000.0f);
    }

    float fovRadians() const { return glm::radians(fovDeg_); }

private:
    glm::vec3 target_   { 0.0f, 0.0f, 0.0f };
    float     distance_ { 0.1f };
    float     yaw_      { glm::radians(35.0f) };
    float     pitch_    { glm::radians(25.0f) };
    float     fovDeg_   { 45.0f };
};

} // namespace waveguide
