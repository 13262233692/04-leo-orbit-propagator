#pragma once
#include <cmath>
#include <array>

namespace leo_propagator {

struct Vector3 {
    double x, y, z;

    Vector3() : x(0.0), y(0.0), z(0.0) {}
    Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vector3 operator+(const Vector3& other) const {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    Vector3 operator*(double scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    Vector3 operator/(double scalar) const {
        return Vector3(x / scalar, y / scalar, z / scalar);
    }

    Vector3& operator+=(const Vector3& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    Vector3& operator*=(double scalar) {
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    double normSquared() const {
        return x*x + y*y + z*z;
    }

    Vector3 normalized() const {
        double n = norm();
        return Vector3(x/n, y/n, z/n);
    }

    double dot(const Vector3& other) const {
        return x*other.x + y*other.y + z*other.z;
    }

    Vector3 cross(const Vector3& other) const {
        return Vector3(
            y*other.z - z*other.y,
            z*other.x - x*other.z,
            x*other.y - y*other.x
        );
    }

    std::array<double, 3> toArray() const {
        return {x, y, z};
    }
};

struct State {
    Vector3 position;
    Vector3 velocity;
    double time;

    State() : position(), velocity(), time(0.0) {}
    State(Vector3 p, Vector3 v, double t) : position(p), velocity(v), time(t) {}
};

}
