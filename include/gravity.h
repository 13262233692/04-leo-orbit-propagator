#pragma once
#include "vector3.h"
#include "constants.h"

namespace leo_propagator {

class GravityModel {
public:
    static Vector3 computeAcceleration(const Vector3& position, bool include_j2 = true);
    static Vector3 computeTwoBodyAcceleration(const Vector3& position);
    static Vector3 computeJ2Acceleration(const Vector3& position);

private:
    static double factorJ2;
};

}
