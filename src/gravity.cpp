#include "gravity.h"
#include <cmath>

namespace leo_propagator {

double factorJ2 = -1.5 * EARTH_J2 * EARTH_MU * EARTH_RADIUS * EARTH_RADIUS;

Vector3 GravityModel::computeAcceleration(const Vector3& position, bool include_j2) {
    Vector3 acc = computeTwoBodyAcceleration(position);
    if (include_j2) {
        acc += computeJ2Acceleration(position);
    }
    return acc;
}

Vector3 GravityModel::computeTwoBodyAcceleration(const Vector3& position) {
    double r = position.norm();
    double r3 = r * r * r;
    double factor = -EARTH_MU / r3;
    return position * factor;
}

Vector3 GravityModel::computeJ2Acceleration(const Vector3& position) {
    double r = position.norm();
    double r2 = r * r;
    double r5 = r2 * r2 * r;
    double z2 = position.z * position.z;
    double ratio = z2 / r2;

    double factor = factorJ2 / r5;

    double ax = factor * position.x * (5.0 * ratio - 1.0);
    double ay = factor * position.y * (5.0 * ratio - 1.0);
    double az = factor * position.z * (5.0 * ratio - 3.0);

    return Vector3(ax, ay, az);
}

}
